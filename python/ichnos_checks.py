"""ICHNOS — high-level validation you run explicitly (not part of a normal
merge): the S=0 leak check, the exogenous-TIP bypass ('Strain 4') check, the
cross-variant shared-value consistency check, and the bypass-model builder
they rely on. These sit ABOVE the core merge (they call build_variant_sbml_string)."""

import libsbml
import tellurium as te

from ichnos_config import (
    TIP_TETR_MODEL, REPORTER_MODEL, VARIANTS,
    CROSS_VARIANT_SHARED_NAMES, _CROSS_VARIANT_ABS_TOL,
)
from ichnos_core import (
    load_model_or_fail, find_species_id_by_name,
    plan_parameter_renames, plan_compartment_renames,
    copy_parameter, copy_species, copy_reaction, copy_rule,
    build_variant_sbml_string,
)
from ichnos_diagnostics import (
    check_unruled_variable_parameters, check_orphan_parameters, find_param_value_anywhere,
)
from ichnos_io import (
    save_merged_sbml, _build_id_to_name_map, _relabel_result_columns, _find_id_by_name,
)

def build_bypass_sbml_string(exogenous_tip_level, save_sbml=False):
    """'Strain 4' logic: bypasses the sensing module ENTIRELY and instead
    holds TIP at a FIXED, externally-supplied level — simulating a
    continuous exogenous TIP supply rather than TIP produced via the
    stress-Hill function. Only merges TIP_TetR_binding + reporter (no
    sensing module at all).

    IMPORTANT (found 2026-08-12): TIP's own degradation lives ONLY inside
    the sensing modules (see Reaction_2 in ERModule/oxidative_module) — the
    TIP_TetR_binding model on its own has no TIP source or sink other than
    binding into TetR_TIP_complex, and complex degradation is a one-way
    sink (TetR_TIP_complex -> null, doesn't release TIP back). So a single
    BOLUS dose (just setting initial concentration) transiently binds and
    then fully drains away with nothing to replenish it — by t=200h any
    dose decays back to the same TIP=0 baseline, making a bolus useless for
    testing steady-state response to exogenous TIP.

    Fix: set the species boundaryCondition=True as well as its
    concentration. This holds TIP CLAMPED at exogenous_tip_level for the
    entire simulation — the binding reaction still consumes it as a
    reactant in the RATE LAW, but SBML boundary species are excluded from
    the ODE integration, so the clamped value never depletes. This properly
    models a continuously-supplied/held external TIP source, giving a real
    non-trivial steady state to compare across doses.
    """
    doc_tetr, m_tetr = load_model_or_fail(TIP_TETR_MODEL)
    doc_reporter, m_reporter = load_model_or_fail(REPORTER_MODEL)

    tip_id = find_species_id_by_name(m_tetr, "TIP")
    tip_species = m_tetr.getSpecies(tip_id)
    tip_species.setInitialConcentration(exogenous_tip_level)
    tip_species.setBoundaryCondition(True)

    reporter_param_renames, reporter_needs_creation = plan_parameter_renames(
        m_tetr, m_reporter, prefix="reporter"
    )
    reporter_compartment_renames = plan_compartment_renames(m_tetr, m_reporter)
    reporter_full_rename_map = {**reporter_param_renames, **reporter_compartment_renames}
    for s in m_reporter.getListOfSpecies():
        copy_species(m_tetr, s, m_reporter)
    for p in m_reporter.getListOfParameters():
        if p.getId() in reporter_needs_creation:
            copy_parameter(m_tetr, p, new_id=reporter_param_renames[p.getId()])
    for i, r in enumerate(m_reporter.getListOfReactions()):
        copy_reaction(m_tetr, r, reporter_full_rename_map, new_id=f"reporter_{i}_{r.getId()}")
    for rule in m_reporter.getListOfRules():
        copy_rule(m_tetr, rule, reporter_full_rename_map)

    check_unruled_variable_parameters(m_tetr, "bypass model")
    check_orphan_parameters(m_tetr, "bypass model")

    sbml_str = libsbml.writeSBMLToString(doc_tetr)
    if save_sbml:
        save_merged_sbml(sbml_str, "bypass", m_tetr, {
            "TIP_TetR_model": TIP_TETR_MODEL, "reporter_module": REPORTER_MODEL,
        })
    return sbml_str


def sanity_check_zero_stress(variant, t_end=200, n_points=200, leak_threshold_frac=0.3):
    """CHECK 1 — 'S = 0' sanity check (POC doc: "απουσία στρες δεν είναι
    επαρκές το TIP"). Runs the variant TWICE:
      (a) with the stress input (S_er / S_ox) forced to 0 — basal/leaky case
      (b) with the stress input forced deep into saturation (1000x its own
          EC50) — the full-activation reference case
    and checks that the basal-case reporter output stays well below the
    fully-activated reference (< leak_threshold_frac, default 30%). This is
    what 'leaky expression shouldn't fully activate the circuit' actually
    means numerically — a plain 'does it run' check wouldn't catch a leaky
    promoter that basically stays half-on all the time.

    Returns a dict with the raw readouts and an overall 'passed' bool. Prints
    a PASS/FAIL summary either way — this is meant to be looked at, not just
    silently trusted.
    """
    stress_name = {"er": "S_er", "ox": "S_ox"}[variant]
    ec50_name = {"er": "EC50_er", "ox": "EC50_ox"}[variant]

    sbml_str = build_variant_sbml_string(variant, save_sbml=False)
    id_to_name = _build_id_to_name_map(sbml_str)
    stress_id = _find_id_by_name(sbml_str, stress_name)
    ec50_id = _find_id_by_name(sbml_str, ec50_name)
    if stress_id is None or ec50_id is None:
        raise ValueError(
            f"Could not find '{stress_name}' and/or '{ec50_name}' by name in the merged "
            f"'{variant}' model — check the sensing module's parameter Names haven't changed."
        )

    def run_at(stress_value):
        r = te.loadSBMLModel(sbml_str)
        r.reset()
        r[stress_id] = stress_value
        result = r.simulate(0, t_end, n_points)
        _relabel_result_columns(result, id_to_name)
        return {name: val for name, val in zip(result.colnames, result[-1])}

    basal = run_at(0.0)
    ec50_value = None
    # grab EC50's own current value to compute a saturating stress level
    r_probe = te.loadSBMLModel(sbml_str)
    ec50_value = r_probe[ec50_id]
    saturating = run_at(ec50_value * 1000.0)

    # 'activation readout' — total reporter signal (sum of all Reporter_* channels)
    reporter_keys = [k for k in basal if "Reporter" in k]
    basal_reporter = sum(basal[k] for k in reporter_keys)
    saturating_reporter = sum(saturating[k] for k in reporter_keys)
    frac = (basal_reporter / saturating_reporter) if saturating_reporter > 0 else float("nan")
    passed = frac < leak_threshold_frac

    print(f"\n=== SANITY CHECK 1 — S=0 (variant={variant}) ===")
    print(f"  Basal (S=0) total reporter signal:        {basal_reporter:.4f}")
    print(f"  Saturating (S=1000×EC50) reporter signal: {saturating_reporter:.4f}")
    print(f"  Basal / saturating ratio: {frac:.3f}  (threshold: < {leak_threshold_frac})")
    print(f"  {'PASS' if passed else 'FAIL'}: basal leak is "
          f"{'well below' if passed else 'NOT below'} the full-activation level.")
    return {
        "variant": variant, "basal": basal, "saturating": saturating,
        "basal_reporter": basal_reporter, "saturating_reporter": saturating_reporter,
        "leak_fraction": frac, "passed": passed,
    }


def sanity_check_exogenous_tip_bypass(tip_levels=(0, 2, 5, 10, 20, 50), t_end=200, n_points=200):
    """CHECK 2 — Exogenous-TIP bypass check ('Strain 4' logic). Runs the
    bypass model (see build_bypass_sbml_string) at several fixed, CLAMPED
    exogenous TIP levels and checks that the feedback loop responds
    MONOTONICALLY — more exogenous TIP should mean more TetR_TIP_complex
    (TIP sequestering TetR away from the promoter), less free TetR_active,
    and therefore MORE total reporter signal (less repression). This
    confirms the downstream TIP→TetR→reporter loop activates correctly
    independent of HOW TIP was produced, not just when it comes from the
    stress-Hill function.

    Prints a table of results and an overall PASS/FAIL on monotonicity (with
    a small tolerance for numerical noise).
    """
    print(f"\n=== SANITY CHECK 2 — Exogenous TIP bypass ('Strain 4' logic) ===")
    rows = []
    for level in tip_levels:
        sbml_str = build_bypass_sbml_string(level, save_sbml=False)
        id_to_name = _build_id_to_name_map(sbml_str)
        r = te.loadSBMLModel(sbml_str)
        r.reset()
        result = r.simulate(0, t_end, n_points)
        _relabel_result_columns(result, id_to_name)
        final = {name: val for name, val in zip(result.colnames, result[-1])}
        reporter_total = sum(v for k, v in final.items() if "Reporter" in k)
        complex_val = next((v for k, v in final.items() if "TetR_TIP_complex" in k), None)
        active_val = next((v for k, v in final.items() if k == "[TetR_active]"), None)
        rows.append({
            "level": level, "reporter_total": reporter_total,
            "TetR_TIP_complex": complex_val, "TetR_active": active_val,
        })
        print(f"  exogenous TIP={level:>6.2f}  ->  TetR_active={active_val:.4f}  "
              f"TetR_TIP_complex={complex_val:.4f}  reporter_total={reporter_total:.4f}")

    reporter_series = [row["reporter_total"] for row in rows]
    tolerance = 1e-6
    monotonic = all(
        reporter_series[i + 1] >= reporter_series[i] - tolerance
        for i in range(len(reporter_series) - 1)
    )
    print(f"  {'PASS' if monotonic else 'FAIL'}: reporter signal is "
          f"{'monotonically non-decreasing' if monotonic else 'NOT monotonic'} with TIP dose.")
    return {"rows": rows, "passed": monotonic}

def sanity_check_qss_speed_robustness(variant, kd_spot_checks=(0.25,), speed_factor=10,
                                       tolerance_frac=0.05, t_end=200, n_points=200):
    """CHECK 3 — QSS speed-robustness spot check. The Option-A sweep (b=4
    fixed, u_w varies -> Kd=u_w/4) only tests the THERMODYNAMIC effect of
    Kd. It does NOT by itself prove that only the ratio u_w/b matters and
    not the absolute speed of binding/unbinding — two models with the same
    Kd but very different absolute rates could in principle behave
    differently if the faster one isn't actually needed to reach QSS, or if
    something non-obvious depends on the raw timescale.

    For each Kd in kd_spot_checks, runs TWO simulations at that same Kd:
      - 'baseline':      b=4  (u_w = Kd*4)               — the normal Option-A point
      - f'{speed_factor}x_faster': b=4*speed_factor  (u_w = Kd*4*speed_factor) —
        SAME Kd, but both rates speed_factor times faster
    and compares final reporter output. If they match within
    tolerance_frac, that CONFIRMS only the ratio (Kd) matters. If they
    DON'T match, that's a real finding: something depends on absolute
    rate, not just Kd.
    """
    print(f"\n=== SANITY CHECK 3 — QSS speed-robustness (variant={variant}, {speed_factor}x spot check) ===")
    all_passed = True
    for kd in kd_spot_checks:
        sbml_str = build_variant_sbml_string(variant, save_sbml=False)
        id_to_name = _build_id_to_name_map(sbml_str)
        doc = libsbml.readSBMLFromString(sbml_str)
        m = doc.getModel()
        b_id = next(p.getId() for p in m.getListOfParameters() if p.getName() == "b")
        # Same name-based lookup as b_id above — SimBiology assigns Kd_TIP_TetR
        # a GUID id on export, keeping "Kd_TIP_TetR" only as the human Name.
        kd_id = next(p.getId() for p in m.getListOfParameters() if p.getName() == "Kd_TIP_TetR")

        def run_at(b_value):
            r = te.loadSBMLModel(sbml_str)
            r.reset()
            r[kd_id] = kd
            r[b_id] = b_value  # rule recomputes u_w = Kd_TIP_TetR * b_value automatically
            result = r.simulate(0, t_end, n_points)
            _relabel_result_columns(result, id_to_name)
            return {name: val for name, val in zip(result.colnames, result[-1])}

        baseline = run_at(4.0)
        faster = run_at(4.0 * speed_factor)

        baseline_reporter = sum(v for k, v in baseline.items() if "Reporter" in k)
        faster_reporter = sum(v for k, v in faster.items() if "Reporter" in k)
        rel_diff = abs(faster_reporter - baseline_reporter) / baseline_reporter if baseline_reporter else float("nan")
        passed = rel_diff < tolerance_frac
        all_passed = all_passed and passed

        print(f"  Kd={kd}: baseline reporter={baseline_reporter:.4f}  "
              f"{speed_factor}x-faster reporter={faster_reporter:.4f}  "
              f"rel.diff={rel_diff:.4f} (threshold {tolerance_frac})")
        print(f"    {'PASS' if passed else 'FAIL'}: {'only the Kd ratio matters, as assumed' if passed else 'absolute rate ALSO matters — Option A alone is insufficient here'}")

    print(f"  OVERALL: {'PASS — QSS assumption holds, Option A sweep is sufficient on its own' if all_passed else 'FAIL — see above, full Option B/C/D parametrization may be needed'}")
    return all_passed

def check_cross_variant_shared_values(variants=None, param_names=None, save_sbml=False):
    """Post-build validation: builds the merged model for EVERY variant and
    verifies that each parameter in CROSS_VARIANT_SHARED_NAMES has the same
    value across all of them.

    WHY THIS EXISTS (2026-08-17): k_deg_TIP is the post-production degradation
    rate of TIP. TIP is produced only inside the sensing modules, so k_deg_TIP
    lives there (Option B — it was removed from the structural binding module,
    which never produced or degraded TIP). But the sensing modules are never
    merged together (ox/er/copper are separate cells → separate builds), so the
    normal in-merge [!] CONFLICT check — which only compares modules that meet
    inside one build — structurally cannot compare k_deg_TIP across variants.
    Yet the value MUST be identical everywhere, because it's the SAME TIP coding
    sequence in all three stresses (confirmed wet lab). This function closes
    exactly that gap: it reads the value from each independently-built variant
    and compares them directly to each other (stricter than the SHARED route,
    which only compares each module against the binding module).

    Raises ValueError on any divergence — this is a hard stop, not a warning,
    because a silent per-variant difference here means the absolute numbers that
    go into the wiki would be inconsistent across stresses for no legitimate
    reason. A value legitimately MISSING from one variant (e.g. that sensing
    module genuinely has no TIP degradation reaction) is also flagged, since
    under the current design every variant is expected to define it.

    Returns a dict {param_name: agreed_value} on success (all variants agree).
    """
    variants = list(variants) if variants is not None else list(VARIANTS)
    param_names = set(param_names) if param_names is not None else set(CROSS_VARIANT_SHARED_NAMES)

    if not param_names:
        return {}

    print(f"\n{'='*60}\nCROSS-VARIANT consistency check\n{'='*60}")
    print(f"  variants checked : {variants}")
    print(f"  params checked   : {sorted(param_names)}")

    # name -> { variant -> (value, scope) }, plus a 'missing' bucket.
    observed = {name: {} for name in param_names}
    missing = {name: [] for name in param_names}

    for v in variants:
        if v not in VARIANTS:
            raise ValueError(f"Unknown variant '{v}', choose from {list(VARIANTS)}")
        # Build (don't re-save) each variant's merged model just to read values.
        # save_sbml defaults to False here so a pure consistency check never
        # rewrites exportsbml\ as a side effect.
        sbml_str = build_variant_sbml_string(v, save_sbml=save_sbml)
        doc = libsbml.readSBMLFromString(sbml_str)
        model = doc.getModel()
        for name in param_names:
            value, scope = find_param_value_anywhere(model, name)
            if value is None:
                missing[name].append(v)
            else:
                observed[name][v] = (value, scope)

    problems = []
    agreed = {}
    for name in sorted(param_names):
        # 1) Missing from one or more variants.
        if missing[name]:
            problems.append(
                f"'{name}' is NOT FOUND in variant(s) {missing[name]} — every "
                f"variant is expected to define it (see CROSS_VARIANT_SHARED_NAMES)."
            )
        # 2) Present but divergent across the variants that do have it.
        values_by_variant = {v: val for v, (val, _sc) in observed[name].items()}
        distinct = sorted(set(values_by_variant.values()))
        # collapse near-equal floats into one representative
        collapsed = []
        for val in distinct:
            if not any(abs(val - c) <= _CROSS_VARIANT_ABS_TOL for c in collapsed):
                collapsed.append(val)
        if len(collapsed) > 1:
            detail = ", ".join(
                f"{v}={val} [{observed[name][v][1]}]"
                for v, val in sorted(values_by_variant.items())
            )
            problems.append(
                f"'{name}' DIVERGES across variants: {detail}. It must be identical "
                f"everywhere (same TIP coding sequence across stresses). Fix the "
                f"differing sensing-module .sbml file(s) so all variants match."
            )
        elif collapsed and not missing[name]:
            agreed[name] = collapsed[0]
            # Show the agreed value + per-variant scope so it's visible where it lives.
            scope_note = "; ".join(
                f"{v}: [{observed[name][v][1]}]" for v in sorted(observed[name])
            )
            print(f"  OK  {name:12s} = {collapsed[0]:<10}  (agrees across {list(observed[name])})")
            print(f"        {scope_note}")

    if problems:
        msg_lines = ["CROSS-VARIANT consistency check FAILED:"]
        for p in problems:
            msg_lines.append(f"  [!] {p}")
        full = "\n".join(msg_lines)
        print("\n" + full)
        print(f"{'='*60}")
        raise ValueError(full)

    print("  → all cross-variant shared parameters agree.")
    print(f"{'='*60}")
    return agreed


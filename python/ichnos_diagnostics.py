"""ICHNOS — low-level diagnostics run DURING a merge: stale-copy detection,
parameter lookup by name, the shared-parameter summary, and the
unruled-variable safety net. These depend only on config (no back-reference
to the merge itself — the higher-level cross-variant check that re-runs the
merge lives in ichnos_checks instead, to keep the dependency graph acyclic)."""

import os
from datetime import datetime

import libsbml

from ichnos_config import SHARED_PARAM_NAMES

def check_for_stale_duplicate(primary_path, alt_path):
    """Two files on disk can share the same basename in different folders:
    'integration\\TIP_TetR_binding.sbml' (the shared, canonical copy the merge
    actually reads) and 'matlab\\tiptetrbinding\\TIP_TetR_binding.sbml' (Ioanna's
    personal working copy). It's expected/normal for her personal copy to be
    ahead while she's mid-edit — this is just a nudge that integration\\ may
    be missing her latest edits, not an error. Best-effort: only compares
    mtimes, so use judgement.
    """
    if not os.path.isfile(primary_path) or not os.path.isfile(alt_path):
        return
    primary_mtime = os.path.getmtime(primary_path)
    alt_mtime = os.path.getmtime(alt_path)
    if alt_mtime > primary_mtime:
        from datetime import datetime
        print(
            f"  [i] NOTE: your personal copy '{alt_path}' was modified more recently "
            f"({datetime.fromtimestamp(alt_mtime)}) than the shared copy this script "
            f"actually uses, '{primary_path}' ({datetime.fromtimestamp(primary_mtime)}). "
            f"If those edits are meant to be in this run, copy them into integration\\ first."
        )


def find_param_value_anywhere(model, target_name):
    """Looks for a parameter by its human-readable NAME (not its raw SBML
    id), both at GLOBAL (model) scope and LOCAL scope (inside each
    reaction's <kineticLaw>). This mirrors find_species_id_by_name's
    name-based approach, for the same reason (see docstring 'GUID
    ιστορία'): tools like COPASI often give global quantities an
    auto-generated SBML id (e.g. something like 'Values[b]' or a numeric/
    GUID-style id) while keeping the human name ('b') as a separate
    attribute — so matching on raw id silently finds nothing even though
    the parameter is right there. Falls back to matching on id itself if a
    parameter has no name set (name attribute is optional in SBML, id is
    not). Returns (value, scope_description) or (None, None) if not found.
    """
    def name_matches(p):
        pname = p.getName()
        return pname == target_name if pname else p.getId() == target_name

    for p in model.getListOfParameters():
        if name_matches(p):
            return p.getValue(), f"GLOBAL (model), id='{p.getId()}'"
    for r in model.getListOfReactions():
        kl = r.getKineticLaw()
        if kl is None:
            continue
        for lp in kl.getListOfParameters():
            if name_matches(lp):
                return lp.getValue(), f"LOCAL to reaction '{r.getId()}', id='{lp.getId()}'"
    return None, None


def print_shared_parameter_summary(model, label):
    """Prints the current value AND scope+id of every watched parameter,
    looked up by NAME (see find_param_value_anywhere), as it ended up in the
    FINAL merged model. A parameter reported as LOCAL is not a problem by
    itself — it just means it's only usable by that one reaction's own
    kinetic law. 'NOT FOUND anywhere' means genuinely missing (wrong file
    loaded, or the name really doesn't exist in this model) — worth
    investigating, unlike a local-scope hit."""
    watch_names = sorted(SHARED_PARAM_NAMES | {
        "b", "u_w", "delta_w", "a_TetR", "k_deg_TetR", "K_R", "n",
        "k_deg_TIP", "P_min",
    })
    print(f"\n--- {label}: merged model parameter check (matched by NAME) ---")
    for name in watch_names:
        value, scope = find_param_value_anywhere(model, name)
        if value is None:
            print(f"  {name:12s} : NOT FOUND anywhere (global or local) by name in the merged model")
        else:
            print(f"  {name:12s} = {value:<10} [{scope}]")
    print("---")


def check_unruled_variable_parameters(model, label):
    """Permanent safety net for the exact bug class copy_rule fixes:
      - every parameter declared constant='false' must have a Rule that
        defines it (a rule-governed quantity like Measured_Ratio_RG).
      - every species declared constant='false' and boundaryCondition='false'
        must be either targeted by a Rule OR be a reactant/product of at
        least one Reaction. A species that's neither has literally nothing
        in the model that ever changes it — it sits frozen at its initial
        value forever, exactly the same silent-freeze failure mode as an
        unruled parameter, just via a different mechanism (missing reaction
        participation instead of missing rule).
    A non-constant value with no defining rule/reaction sits at whatever
    initial value it was given, frozen forever — no error, no dangling-
    symbol warning, simulation runs fine, just silently wrong. Prints a loud
    warning for each offender found; does not raise, since a deliberately
    free variable is technically valid SBML too — but in practice this is
    the bug, not the exception, so it's worth seeing every time.
    """
    ruled_variables = {
        rule.getVariable() for rule in model.getListOfRules() if hasattr(rule, "getVariable")
    }
    reacted_species_ids = set()
    for r in model.getListOfReactions():
        for ref in r.getListOfReactants():
            reacted_species_ids.add(ref.getSpecies())
        for ref in r.getListOfProducts():
            reacted_species_ids.add(ref.getSpecies())

    offenders = []
    for p in model.getListOfParameters():
        if not p.getConstant() and p.getId() not in ruled_variables:
            offenders.append(f"parameter '{p.getName() or p.getId()}' (id={p.getId()}) — no Rule defines it")
    for s in model.getListOfSpecies():
        if s.getConstant() or s.getBoundaryCondition():
            continue
        if s.getId() not in ruled_variables and s.getId() not in reacted_species_ids:
            offenders.append(
                f"species '{s.getName() or s.getId()}' (id={s.getId()}) — no Rule defines it "
                f"AND it's not a reactant/product of any Reaction"
            )
    if offenders:
        print(f"\n  [!] UNRULED VARIABLE(S) in {label}: the following are declared "
              f"constant='false' (meant to be computed) but nothing in the model updates them — "
              f"they will sit frozen at their initial value for the entire simulation:")
        for o in offenders:
            print(f"      - {o}")
    return offenders


def check_missing_units(model, label):
    """Reports every parameter that has NO units attribute set. In a model
    where concentrations are nM and rates are 1/h, a unitless parameter like
    EC50=271 is ambiguous (nM? uM?) — a documentation/interpretation hazard,
    not a simulation error (roadrunner ignores units at runtime). Target: 0.
    Returns the list of offending parameter ids."""
    missing = [p.getId() for p in model.getListOfParameters() if not p.isSetUnits()]
    if missing:
        print(f"\n  [!] {len(missing)} parameter(s) WITHOUT units in {label}: {missing}")
    else:
        print(f"\n  [i] units check ({label}): all parameters have units declared.")
    return missing

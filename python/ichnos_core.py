"""ICHNOS — the merge core. Reads the separate SBML submodels and combines
them IN MEMORY (name-based matching, never id-based — see GUID history in the
config docstring). This is the part that actually produces the merged model;
everything else supports or inspects it."""

import os

import libsbml

from ichnos_config import (
    TIP_TETR_MODEL, REPORTER_MODEL, TIP_TETR_MODEL_ALT_CHECK, VARIANTS,
    SHARED_PARAM_NAMES, SHARED_PARAM_TOLERANCE,
)
from ichnos_diagnostics import (
    check_for_stale_duplicate, check_unruled_variable_parameters,
    print_shared_parameter_summary,
)
from ichnos_io import save_merged_sbml

def load_model_or_fail(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"Δεν βρέθηκε το αρχείο: {path}\n"
            f"  (working directory: {os.getcwd()})\n"
            f"  Έλεγξε ότι το path είναι σωστό/απόλυτο, ή ότι τρέχεις το script απ' τον σωστό φάκελο."
        )
    doc = libsbml.readSBMLFromFile(path)
    model = doc.getModel()
    if model is None:
        errs = "\n".join(doc.getError(i).getMessage() for i in range(doc.getNumErrors()))
        raise ValueError(f"Το αρχείο βρέθηκε αλλά δεν parse-άρισε σωστά σε SBML model: {path}\n{errs}")
    return doc, model


def find_species_id_by_name(model, name):
    for s in model.getListOfSpecies():
        if s.getName() == name:
            return s.getId()
    raise ValueError(f"No species named '{name}' found in model {model.getId()}")


def find_param_id_by_name(model, name):
    for p in model.getListOfParameters():
        if p.getName() == name:
            return p.getId()
    return None


def _param_display_name(p):
    """The meaningful identity of a parameter is its Name (e.g. 'mu', 'n'),
    NOT its raw SBML id — confirmed 2026-08-12 that this COPASI-exported
    model gives every parameter a GUID-style id (e.g.
    'mwfc0b18ce_ba8a_4f91_bb9e_838d6f352026') and keeps the human name as a
    separate attribute. Falls back to id only if Name is genuinely unset."""
    return p.getName() or p.getId()


def rename_in_ast(node, mapping):
    if node is None:
        return
    if node.isName() and node.getName() in mapping:
        node.setName(mapping[node.getName()])
    for i in range(node.getNumChildren()):
        rename_in_ast(node.getChild(i), mapping)


def resolve_compartment_id(dest_model, source_model, source_compartment_id):
    """Given a compartment id as used by a species/reaction in source_model,
    returns the compartment id to actually use once merged into dest_model.

    Same class of problem as species/parameter ids (see 'GUID ιστορία' /
    _param_display_name): a compartment called 'cell' in the reporter module
    might not exist under that exact id in Ioanna's model, even if
    conceptually it's the same physical compartment. Resolution order:
      1. dest already has a compartment with this exact id -> reuse it.
      2. dest has a DIFFERENT-id compartment with the SAME NAME -> reuse that
         one's id instead (this is the expected fix for the
         'references unknown compartment cell' error).
      3. Neither matches -> copy the compartment definition from source into
         dest, preserving its id, but WARN loudly — this likely means the two
         files model compartments differently (e.g. single-compartment vs
         multi-compartment) and deserves a human look, not a silent patch.
    """
    if dest_model.getCompartment(source_compartment_id) is not None:
        return source_compartment_id

    src_c = source_model.getCompartment(source_compartment_id)
    src_name = src_c.getName() if src_c is not None else None
    if src_name:
        for c in dest_model.getListOfCompartments():
            if c.getName() == src_name:
                return c.getId()

    if src_c is None:
        raise ValueError(
            f"Compartment '{source_compartment_id}' is referenced by a species/reaction but "
            f"is not even defined in its own source model — the source .sbml file itself is "
            f"malformed."
        )
    print(
        f"  [!] COMPARTMENT MISSING: '{source_compartment_id}'"
        + (f" (name='{src_name}')" if src_name else "")
        + f" not found in the destination model by id OR by name — creating it by copying "
        f"the definition from the source file as-is. This likely means the two .sbml files "
        f"model compartments differently (e.g. different names for what should be the same "
        f"physical compartment, or a genuinely separate compartment) — worth a human check, "
        f"not just trusting this fallback."
    )
    new_c = dest_model.createCompartment()
    new_c.setId(src_c.getId())
    new_c.setName(src_c.getName())
    new_c.setConstant(src_c.getConstant())
    if src_c.isSetSize():
        new_c.setSize(src_c.getSize())
    return src_c.getId()


def plan_compartment_renames(dest_model, source_model):
    """Resolves EVERY compartment defined in source_model against dest_model
    (via resolve_compartment_id) and returns {old_id: resolved_id}.

    This must be merged into the SAME id_rename_map that's passed to
    copy_reaction — not just used for species. Reason: resolve_compartment_id
    already fixes the compartment attribute on copied SPECIES correctly, but
    a reaction's kinetic-law MATH can also reference a compartment id
    directly (e.g. a volume-scaling term like '.../cell', common in
    COPASI-exported concentration rate laws). Without this rename applied to
    the reaction math too, the kinetic law keeps the literal source
    compartment id, which may not exist at all in the merged model (only its
    resolved/GUID counterpart does) — roadrunner then fails at load time with
    something like \"symbol 'cell' is not physically stored... it either does
    not exist or is defined by an assignment rule\", which is a confusing
    error for what is actually just a missed rename.
    """
    return {
        c.getId(): resolve_compartment_id(dest_model, source_model, c.getId())
        for c in source_model.getListOfCompartments()
    }


def plan_parameter_renames(dest_model, source_model, prefix, skip_names=frozenset()):
    """Decides what to do with every parameter in source_model, matched by
    NAME against dest_model's CURRENT global parameters (not raw id — see
    _param_display_name). Three outcomes per parameter:

    - name in skip_names: not copied at all, and left OUT of the rename map
      entirely (so any reference to it in the source's own reactions is left
      as the literal original name). NOT currently used by any call site as
      of 2026-08-13 — 'P' used to be handled this way (relying on a literal
      id match that turned out to never hold, since every id in these
      COPASI-exported files is a GUID) but is now correctly handled via
      SHARED_PARAM_NAMES below instead, same as 'mu'. skip_names is kept as
      a general mechanism for a genuinely different future need (something
      that should be dropped entirely, not unified) — it is NOT how 'P' is
      currently resolved, despite what an older version of this docstring
      said.
    - name in SHARED_PARAM_NAMES (e.g. 'mu', 'P'): if dest already has a global
      parameter with this name (regardless of its — possibly GUID — id),
      this source parameter is the SAME global quantity. It is NOT copied as
      a second parameter; instead its id is mapped onto dest's existing id,
      so its own reactions correctly point at the one shared copy. A value
      mismatch is reported (not silently dropped). If dest does NOT yet have
      it, this source's copy becomes the first (and future submodules in the
      same merge will then match against it).
    - anything else: if dest already has a DIFFERENT global parameter with
      the same name (a genuine, likely accidental, collision — e.g. two
      unrelated Hill coefficients both named 'n'), this copy gets a fresh,
      guaranteed-unique id and a loud warning. Otherwise it keeps its own
      original id unchanged (GUID ids essentially never collide raw, but we
      still register it in the name index so subsequent submodules in the
      same merge are checked against it too).

    Returns (rename_map, needs_creation):
      rename_map:     {old_id: final_id} — feed the FULL merged dict (this +
                       any species rename) into copy_reaction's id_rename_map
                       so kinetic-law math ends up pointing at the right ids.
      needs_creation:  set of old_ids that should actually be passed to
                       copy_parameter() to create a new SBML parameter.
                       old_ids NOT in this set (the 'already shared, matched
                       by name' case) must NOT be created — creating them
                       would just silently duplicate an existing global
                       parameter under a second id.
    """
    dest_name_index = {}
    for p in dest_model.getListOfParameters():
        dest_name_index[_param_display_name(p)] = p.getId()

    rename_map = {}
    needs_creation = set()
    for p in source_model.getListOfParameters():
        old_id = p.getId()
        pname = _param_display_name(p)
        if pname in skip_names:
            continue
        if pname in SHARED_PARAM_NAMES:
            if pname in dest_name_index:
                existing_id = dest_name_index[pname]
                rename_map[old_id] = existing_id
                existing = dest_model.getParameter(existing_id)
                if existing is not None and abs(existing.getValue() - p.getValue()) > SHARED_PARAM_TOLERANCE:
                    print(
                        f"  [!] CONFLICT: shared parameter '{pname}' already = {existing.getValue()} "
                        f"in the destination model (id='{existing_id}'), but {prefix} defines it as "
                        f"{p.getValue()} (id='{old_id}'). Keeping the destination's value "
                        f"({existing.getValue()}); the {prefix} value is being DISCARDED. "
                        f"Fix this in the .sbml files if that's not intended."
                    )
            else:
                rename_map[old_id] = old_id
                needs_creation.add(old_id)
                dest_name_index[pname] = old_id
            continue
        if pname in dest_name_index:
            new_id = f"{prefix}_{old_id}"
            while dest_model.getParameter(new_id) is not None:
                new_id = f"{new_id}_"
            print(
                f"  [!] NAME COLLISION: parameter named '{pname}' from {prefix} (id='{old_id}') "
                f"has the SAME NAME as an already-merged parameter (likely a DIFFERENT quantity, "
                f"e.g. two different Hill coefficients both called 'n') — renaming this copy to "
                f"id='{new_id}' so its value/meaning isn't silently confused with the existing one. "
                f"Its reactions are renamed to match automatically."
            )
            rename_map[old_id] = new_id
            needs_creation.add(old_id)
        else:
            rename_map[old_id] = old_id
            needs_creation.add(old_id)
            dest_name_index[pname] = old_id
    return rename_map, needs_creation


def copy_parameter(dest_model, param, new_id):
    """Creates param in dest_model under new_id. Caller (via
    plan_parameter_renames's needs_creation set) is responsible for only
    calling this when new_id is NOT already taken — this function no longer
    silently no-ops on a collision, since with name-based planning upstream,
    reaching an actual id collision here means something skipped the plan
    and is a real bug worth a loud failure rather than silently discarding
    data."""
    if dest_model.getParameter(new_id) is not None:
        raise AssertionError(
            f"copy_parameter called with new_id='{new_id}' which already exists in the "
            f"destination model. This should be impossible if plan_parameter_renames was used "
            f"correctly — only call copy_parameter for ids in its 'needs_creation' set."
        )
    p = dest_model.createParameter()
    p.setId(new_id)
    p.setName(param.getName())
    p.setValue(param.getValue())
    p.setConstant(param.getConstant())


def copy_species(dest_model, species, source_model, new_id=None):
    pid = new_id or species.getId()
    if dest_model.getSpecies(pid) is not None:
        # Same class of bug as the parameter id collision: creating a second
        # species with an id that already exists produces invalid/ambiguous
        # SBML. Currently only the reporter module's species are copied this
        # way, and only after the TIP id has already been unified — so this
        # should never fire in practice, but fail loudly instead of emitting
        # broken SBML if a future submodel change introduces a clash.
        raise ValueError(
            f"Species id collision: '{pid}' already exists in the destination model. "
            f"Rename it in the source .sbml file, or extend copy_species with the same "
            f"rename-planning approach used for parameters (plan_parameter_renames)."
        )
    resolved_compartment = resolve_compartment_id(dest_model, source_model, species.getCompartment())
    s = dest_model.createSpecies()
    s.setId(new_id or species.getId())
    s.setName(species.getName())
    s.setCompartment(resolved_compartment)
    s.setInitialConcentration(species.getInitialConcentration())
    s.setConstant(species.getConstant())
    s.setBoundaryCondition(species.getBoundaryCondition())
    s.setHasOnlySubstanceUnits(species.getHasOnlySubstanceUnits())


def copy_reaction(dest_model, reaction, id_rename_map, new_id):
    r = dest_model.createReaction()
    r.setId(new_id)
    r.setReversible(reaction.getReversible())
    r.setFast(False)
    for i in range(reaction.getNumReactants()):
        ref = reaction.getReactant(i)
        nref = r.createReactant()
        nref.setSpecies(id_rename_map.get(ref.getSpecies(), ref.getSpecies()))
        nref.setStoichiometry(ref.getStoichiometry())
        nref.setConstant(True)
    for i in range(reaction.getNumProducts()):
        ref = reaction.getProduct(i)
        nref = r.createProduct()
        nref.setSpecies(id_rename_map.get(ref.getSpecies(), ref.getSpecies()))
        nref.setStoichiometry(ref.getStoichiometry())
        nref.setConstant(True)
    kl_src = reaction.getKineticLaw()
    if kl_src is not None:
        math_copy = kl_src.getMath().deepCopy()
        # Single pass over the FULL mapping — rename_in_ast already walks the
        # whole AST and substitutes every key it finds. Looping and calling it
        # once per mapping entry (as before) was redundant when the mapping
        # only ever had 1 entry, and would be actively wrong for a mapping
        # with 2+ entries that chain (e.g. A->B and B->C would cascade into
        # A->C on a second pass). One call is correct and sufficient.
        rename_in_ast(math_copy, id_rename_map)
        kl = r.createKineticLaw()
        kl.setMath(math_copy)


def copy_rule(dest_model, rule, id_rename_map):
    """Copies an SBML Rule (AssignmentRule / RateRule / AlgebraicRule) into
    dest_model, renaming BOTH its target variable and its math via
    id_rename_map — same treatment as copy_reaction's kinetic law.

    ADDED 2026-08-13: this was MISSING entirely. The original merge only
    copied species/parameters/reactions (three loops, three copy_* helpers)
    — never rules. A parameter with constant='false' whose value is meant to
    be computed by a rule (e.g. reporter's Measured_Ratio_RG, Ratio_RG_FRET,
    Observed_Green, Total_red_pool, b_fret) got copied as a plain static
    parameter (copy_parameter doesn't care WHY something is non-constant,
    it just copies id/name/value/constant-flag) — its INITIAL value, frozen
    forever, since nothing ever recomputes it. No error, no dangling-symbol
    warning (the parameter genuinely exists) — the simulation runs and looks
    fine, it's just silently wrong for anything reading that parameter.

    The variable's OWN id-renaming was already handled correctly by
    plan_parameter_renames/copy_parameter (rule targets are ordinary
    parameters as far as that logic is concerned) — this only adds the
    missing rule itself, pointed at the same (possibly renamed) variable id.
    """
    kind = rule.getElementName()  # 'assignmentRule', 'rateRule', 'algebraicRule'
    variable = rule.getVariable() if hasattr(rule, "getVariable") else None
    new_variable = id_rename_map.get(variable, variable) if variable else None

    # SBML forbids more than one Rule targeting the same variable (L2V4 §4.11).
    # This can happen here specifically for SHARED_PARAM_NAMES variables
    # (mu, P): if a submodule defines ITS OWN rule for what gets recognized
    # as 'the same' shared variable, that rule would collide with dest's
    # existing one for that variable. Keeping dest's existing rule (i.e. NOT
    # copying this one) is almost certainly correct — shared params exist
    # precisely to be unified — but it must be loud, not silent.
    if new_variable is not None:
        existing_rules = [
            r for r in dest_model.getListOfRules() if r.getVariable() == new_variable
        ]
        if existing_rules:
            print(
                f"  [!] DUPLICATE RULE: variable '{new_variable}' already has a rule in the "
                f"destination model — this incoming rule is NOT copied (SBML forbids two rules "
                f"on the same variable). If this is a shared variable (e.g. mu, P), keeping "
                f"dest's existing rule is expected; if not, the two definitions genuinely "
                f"conflict and need a look."
            )
            return None

    math_copy = rule.getMath().deepCopy()
    rename_in_ast(math_copy, id_rename_map)

    if rule.isAssignment():
        new_rule = dest_model.createAssignmentRule()
    elif rule.isRate():
        new_rule = dest_model.createRateRule()
    elif rule.isAlgebraic():
        new_rule = dest_model.createAlgebraicRule()
    else:
        raise ValueError(f"Unsupported rule type '{kind}' — extend copy_rule to handle it.")

    if new_variable is not None:
        new_rule.setVariable(new_variable)
    new_rule.setMath(math_copy)
    return new_rule


def build_variant_sbml_string(variant, save_sbml=True):
    """Returns an SBML string for the requested variant, built fresh from the
    separate source files. By default ALSO writes it to exportsbml\\ (a
    'latest' copy plus a timestamped, never-overwritten archive copy with a
    manifest — see save_merged_sbml) for downstream use in sensitivity
    analysis / the comparison circuit. Pass save_sbml=False to skip that and
    keep the old in-memory-only behavior."""
    cfg = VARIANTS[variant]

    check_for_stale_duplicate(TIP_TETR_MODEL, TIP_TETR_MODEL_ALT_CHECK)

    doc_tetr, m_tetr = load_model_or_fail(TIP_TETR_MODEL)
    doc_sensing, m_sensing = load_model_or_fail(cfg["sensing_file"])
    doc_reporter, m_reporter = load_model_or_fail(REPORTER_MODEL)

    tip_id = find_species_id_by_name(m_tetr, "TIP")
    sensing_tip_id = find_species_id_by_name(m_sensing, cfg["tip_name"])

    # --- sensing module reactions/params merge straight into Ioanna's model ---
    # Plan collisions (matched by NAME, not raw id — see plan_parameter_renames)
    # BEFORE copying anything, so both the parameters themselves AND the
    # reaction math that references them get renamed consistently. E.g. if
    # the sensing module also has a parameter NAMED 'n' (its own Hill
    # coefficient, e.g. 1.7 for ox / 2.13 for er) that is a DIFFERENT
    # quantity from Ioanna's 'n' (TetR cooperativity, =4), this renames the
    # sensing one to a fresh id — instead of silently discarding it and
    # having the sensing reaction accidentally use Ioanna's n=4. Conversely
    # if the sensing module has its own 'mu' (SHARED_PARAM_NAMES), it's
    # recognized as the SAME global mu and pointed at Ioanna's existing one
    # instead of creating a duplicate.
    sensing_param_renames, sensing_needs_creation = plan_parameter_renames(
        m_tetr, m_sensing, prefix=f"sensing_{variant}"
    )
    # Also resolve any compartment ids the sensing module's own reactions
    # might reference directly in their kinetic-law math (volume-scaling
    # terms etc.) — even though we don't copy the sensing module's species,
    # its reaction MATH can still mention a compartment id by itself. See
    # plan_compartment_renames docstring for why this is needed alongside
    # (not instead of) the species-level resolve_compartment_id calls.
    sensing_compartment_renames = plan_compartment_renames(m_tetr, m_sensing)
    sensing_full_rename_map = {
        sensing_tip_id: tip_id,
        **sensing_param_renames,
        **sensing_compartment_renames,
    }

    for p in m_sensing.getListOfParameters():
        old_id = p.getId()
        if old_id in sensing_needs_creation:
            copy_parameter(m_tetr, p, new_id=sensing_param_renames[old_id])
    for i, r in enumerate(m_sensing.getListOfReactions()):
        copy_reaction(m_tetr, r, sensing_full_rename_map, new_id=f"sensing_{i}_{r.getId()}")
    for rule in m_sensing.getListOfRules():
        copy_rule(m_tetr, rule, sensing_full_rename_map)

    # --- reporter module merge, skipping its own local P (use Ioanna's dynamic P instead) ---
    # 'P' used to be excluded via skip_names on the (wrong) assumption that
    # leaving it un-renamed would make it fall through to Ioanna's real P by
    # a literal id match ("P" == "P"). Now that we know every id in this
    # model is a GUID, that assumption never held — 'P' is handled instead
    # via SHARED_PARAM_NAMES below (same mechanism as 'mu'): resolved BY
    # NAME to Ioanna's actual (rule-governed) P parameter id.
    reporter_param_renames, reporter_needs_creation = plan_parameter_renames(
        m_tetr, m_reporter, prefix="reporter"
    )
    # Resolve reporter's compartment(s) too, and fold into the SAME rename
    # map used for its reactions — this is the fix for the 'cell' load error:
    # copy_species below correctly resolves each species' compartment
    # attribute, but that alone doesn't touch a compartment id mentioned
    # directly inside a reaction's kinetic-law math (e.g. a volume-scaling
    # term) — this does.
    reporter_compartment_renames = plan_compartment_renames(m_tetr, m_reporter)
    reporter_full_rename_map = {**reporter_param_renames, **reporter_compartment_renames}

    for s in m_reporter.getListOfSpecies():
        copy_species(m_tetr, s, m_reporter)
    for p in m_reporter.getListOfParameters():
        old_id = p.getId()
        if old_id in reporter_needs_creation:
            copy_parameter(m_tetr, p, new_id=reporter_param_renames[old_id])
    for i, r in enumerate(m_reporter.getListOfReactions()):
        copy_reaction(m_tetr, r, reporter_full_rename_map, new_id=f"reporter_{i}_{r.getId()}")
    for rule in m_reporter.getListOfRules():
        copy_rule(m_tetr, rule, reporter_full_rename_map)

    check_unruled_variable_parameters(m_tetr, f"variant={variant}")
    print_shared_parameter_summary(m_tetr, f"variant={variant}")

    sbml_str = libsbml.writeSBMLToString(doc_tetr)
    if save_sbml:
        source_paths = {
            "TIP_TetR_model": TIP_TETR_MODEL,
            "sensing_module": cfg["sensing_file"],
            "reporter_module": REPORTER_MODEL,
        }
        try:
            save_merged_sbml(sbml_str, variant, m_tetr, source_paths)
        except OSError as e:
            # Saving the merged SBML to disk is a convenience for
            # sensitivity analysis / the comparison circuit, NOT required for
            # the simulation itself — never let a filesystem hiccup (OneDrive
            # sync locks, permissions, etc.) abort the actual run. Warn and
            # keep going with the in-memory model.
            print(
                f"  [!] WARNING: could not save merged SBML to disk ({e}). "
                f"Continuing with the in-memory model anyway — nothing archived "
                f"for variant='{variant}' this run. Re-run once the folder is "
                f"accessible if you need the archive/manifest."
            )
    return sbml_str


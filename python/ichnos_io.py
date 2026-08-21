"""ICHNOS — file I/O: saving merged SBML (latest + archive + manifest),
and display helpers that turn GUID ids into readable Names for plots/prints."""

import os
import json
from datetime import datetime

import libsbml

from ichnos_config import EXPORT_SBML_DIR, SHARED_PARAM_NAMES
from ichnos_diagnostics import find_param_value_anywhere

def _ensure_dir(path):
    """Best-effort directory creation. On this machine, os.path.isdir() has
    proven unreliable for a folder inside OneDrive (Files-On-Demand can
    present it as a cloud placeholder / reparse point that Python's stat
    calls don't always resolve consistently) — so we no longer trust isdir()
    to decide anything. Just attempt makedirs and swallow FileExistsError
    unconditionally: if something is there, proceed; if it turns out not to
    be a real usable directory, the subsequent open(..., 'w') call will fail
    with its own clear error instead."""
    try:
        os.makedirs(path, exist_ok=True)
    except FileExistsError:
        pass


def _collect_shared_param_snapshot(model):
    watch_names = sorted(SHARED_PARAM_NAMES | {
        "b", "u_w", "delta_w", "a_TetR", "k_deg_TetR", "K_R", "n",
        "k_deg_TIP", "P_min",
    })
    snapshot = {}
    for name in watch_names:
        value, scope = find_param_value_anywhere(model, name)
        snapshot[name] = {"value": value, "scope": scope}
    return snapshot


def save_merged_sbml(sbml_str, variant, m_tetr, source_paths):
    """Writes the merged SBML to two places:
      - exportsbml\\merged_<variant>.sbml        — 'latest' pointer, overwritten every run
      - exportsbml\\archive\\merged_<variant>_<timestamp>.sbml — never overwritten, plus a
        sidecar .json manifest recording exactly which source .sbml files (and their
        mtimes) and which key parameter values went into THIS specific merge.

    The archive copy exists so that downstream work (sensitivity analysis,
    the no-feedback comparison circuit, etc.) can always point back at the
    exact merged model a given result came from, instead of relying on
    whatever 'merged_<variant>.sbml' happens to contain right now.
    Returns (latest_path, archive_path, manifest_path).
    """
    import json
    from datetime import datetime

    _ensure_dir(EXPORT_SBML_DIR)
    archive_dir = os.path.join(EXPORT_SBML_DIR, "archive")
    _ensure_dir(archive_dir)

    latest_path = os.path.join(EXPORT_SBML_DIR, f"merged_{variant}.sbml")
    with open(latest_path, "w", encoding="utf-8") as f:
        f.write(sbml_str)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    archive_path = os.path.join(archive_dir, f"merged_{variant}_{timestamp}.sbml")
    with open(archive_path, "w", encoding="utf-8") as f:
        f.write(sbml_str)

    manifest = {
        "variant": variant,
        "built_at": datetime.now().isoformat(timespec="seconds"),
        "source_files": {
            label: {
                "path": path,
                "modified": datetime.fromtimestamp(os.path.getmtime(path)).isoformat(timespec="seconds")
                if os.path.isfile(path) else None,
            }
            for label, path in source_paths.items()
        },
        "key_parameters_in_merged_model": _collect_shared_param_snapshot(m_tetr),
    }
    manifest_path = archive_path.replace(".sbml", ".manifest.json")
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)

    print(f"  Merged SBML (latest):  {latest_path}")
    print(f"  Merged SBML (archive): {archive_path}")
    print(f"  Manifest:              {manifest_path}")
    return latest_path, archive_path, manifest_path


def _build_id_to_name_map(sbml_str):
    """Maps every species/parameter SBML id to its human-readable Name (see
    _param_display_name — same GUID-vs-Name situation applies here). Used
    purely for display (plot legends, printed values) — never for the merge
    logic itself, which must keep using real ids."""
    doc = libsbml.readSBMLFromString(sbml_str)
    m = doc.getModel()
    id_to_name = {}
    for s in m.getListOfSpecies():
        id_to_name[s.getId()] = s.getName() or s.getId()
    for p in m.getListOfParameters():
        id_to_name[p.getId()] = p.getName() or p.getId()
    return id_to_name


def _relabel_result_columns(result, id_to_name):
    """Rewrites a roadrunner simulate() result's column headers from raw
    SBML ids (e.g. '[mw7303685c_...]') to their human-readable Names (e.g.
    '[TetR_active]'), so plot legends and printed 'Final values' are
    actually readable. roadrunner's NamedArray.colnames is directly
    settable — this doesn't touch the underlying data, only the labels."""
    new_colnames = []
    for col in result.colnames:
        if col == "time":
            new_colnames.append(col)
            continue
        inner = col.strip("[]")
        new_colnames.append(f"[{id_to_name.get(inner, inner)}]")
    result.colnames = new_colnames
    return result


def _find_id_by_name(sbml_str, target_name):
    """Reverse of _build_id_to_name_map: given a human-readable Name, finds
    the actual SBML id to use with roadrunner's r[id] getter/setter (which
    operates on real ids, not Names)."""
    doc = libsbml.readSBMLFromString(sbml_str)
    m = doc.getModel()
    for p in m.getListOfParameters():
        if (p.getName() or p.getId()) == target_name:
            return p.getId()
    for s in m.getListOfSpecies():
        if (s.getName() or s.getId()) == target_name:
            return s.getId()
    return None


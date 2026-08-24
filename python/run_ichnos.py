"""ICHNOS — entry point. Run the merge + simulation for one or both stress
variants (er/ox), with optional sanity checks and export.

Usage:
    python run_ichnos.py            # both variants (er, ox), saves exportsbml/
    python run_ichnos.py er         # only ER-stress variant
    python run_ichnos.py ox         # only oxidative-stress variant
    python run_ichnos.py --no-save  # don't write anything to exportsbml/
    python run_ichnos.py --sanity   # run the S=0 and exogenous-TIP sanity checks

This file is deliberately thin: the real work lives in ichnos_core (merge),
ichnos_diagnostics (in-merge checks), ichnos_io (save/labels), and
ichnos_checks (explicit validation). See ichnos_config for paths/constants.
"""

import sys

from ichnos_config import VARIANTS
from ichnos_core import build_variant_sbml_string
from ichnos_io import _build_id_to_name_map, _relabel_result_columns
from ichnos_checks import (
    sanity_check_zero_stress, sanity_check_exogenous_tip_bypass,
    sanity_check_qss_speed_robustness, check_cross_variant_shared_values,
)

import libsbml
import tellurium as te


def run_variant(variant, t_end=200, n_points=500, plot=True, save_sbml=True):
    print(f"\n{'='*60}\nRunning variant='{variant}'\n{'='*60}")
    sbml_str = build_variant_sbml_string(variant, save_sbml=save_sbml)
    id_to_name = _build_id_to_name_map(sbml_str)
    r = te.loadSBMLModel(sbml_str)

    r.reset()
    result = r.simulate(0, t_end, n_points)
    _relabel_result_columns(result, id_to_name)

    print("Final values:")
    for name, val in zip(result.colnames, result[-1]):
        print(f"  {name}: {val:.4f}")

    if plot:
        r.plot(result, title=f"ICHNOS — variant={variant}")

    return r, result


if __name__ == "__main__":
    args = sys.argv[1:]
    save_sbml = "--no-save" not in args
    run_sanity = "--sanity" in args
    args = [a for a in args if a not in ("--no-save", "--sanity")]
    requested = args or ["er", "ox"]

    if run_sanity:
            all_passed = True
    for v in requested:
            if v not in VARIANTS:
                print(f"Unknown variant '{v}', choose from {list(VARIANTS)}")
                sys.exit(1)
            result1 = sanity_check_zero_stress(v)
            all_passed = all_passed and result1["passed"]
            qss_passed = sanity_check_qss_speed_robustness(v, kd_spot_checks=(0.01, 0.25, 5.0), speed_factor=10)
            all_passed = all_passed and qss_passed
    result2 = sanity_check_exogenous_tip_bypass()
    all_passed = all_passed and result2["passed"]
    print(f"\n{'='*60}\nOVERALL: {'ALL CHECKS PASSED' if all_passed else 'SOME CHECKS FAILED — see above'}\n{'='*60}")
    sys.exit(0 if all_passed else 1)
    

    for v in requested:
        if v not in VARIANTS:
            print(f"Unknown variant '{v}', choose from {list(VARIANTS)} (or add --no-save to skip writing exportsbml\\*.sbml)")
            sys.exit(1)

    # Cross-variant consistency gate (see check_cross_variant_shared_values):
    # only meaningful when >1 variant is in play — a single-variant run has
    # nothing to compare against. Runs BEFORE any simulation so an inconsistent
    # k_deg_TIP across stresses stops the run instead of quietly producing
    # per-stress-inconsistent numbers. Uses save_sbml=False so this pure check
    # never rewrites exportsbml\ (the real per-variant save happens in
    # run_variant below).
    if len(requested) > 1:
        try:
            check_cross_variant_shared_values(variants=requested, save_sbml=False)
        except ValueError:
            print("\n[!] Aborting before simulation due to cross-variant inconsistency above.")
            sys.exit(1)
    else:
        print(f"\n[i] Single variant ('{requested[0]}') requested — skipping cross-variant "
              f"consistency check (nothing to compare against). Run with ≥2 variants to check.")

    for v in requested:
        run_variant(v, save_sbml=save_sbml)

        
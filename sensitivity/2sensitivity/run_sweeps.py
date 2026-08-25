"""
ICHNOS — sensitivity sweeps (Φάση 3)

Τρέχει τα τρία OAT sweeps (K_R, Kd_TIP-TetR μέσω u_w, n) πάνω στα merged
μοντέλα, για κάθε variant, και υπολογίζει τα τρία μετρικά ανά καμπύλη.

    python run_sweeps.py

Έξοδος: sweep_results.csv  +  εκτύπωση έτοιμη για αντιγραφή στο workbook.

ΣΗΜΕΙΩΣΕΙΣ
- Έξοδος του μοντέλου: Observed_Green (ΟΧΙ Measured_Ratio_RG — ο tandem timer
  είναι εξ ορισμού ανεξάρτητος του ρυθμού παραγωγής, Khmelinskii Note 4).
- Κάθε γραμμή αποτελέσματος = ΜΙΑ ολόκληρη καμπύλη dose-response (8 τρεξίματα).
- Το S ορίζεται προγραμματιστικά· στα .sbml παραμένει 0.
"""

import csv
import numpy as np
import roadrunner
from scipy.optimize import curve_fit

MERGED = {"ox": "exportsbml/merged_ox.sbml",
          "er": "exportsbml/merged_er.sbml"}

# ids επαληθευμένα στα merged αρχεία (2026-08-22).
# Τα GUID ids ανήκουν στο μοντέλο της Ιωάννας, τα απλά στα sensing/reporter modules.
IDS = {
    "ox": {"S": "S_ox", "out": "Observed_Green"},
    "er": {"S": "S_er", "out": "Observed_Green"},
}
P_KR  = "mw6cc4754e_4133_4feb_a065_4df7424c2bb5"   # K_R   = 0.44 nM
P_UW  = "mwcd6e7518_be87_4a94_9eb4_d15c7b443e24"   # u_w   = 1 1/h   (Kd = u_w/4)
P_N   = "mwbe341ad3_9a73_4ff9_aa05_5c9969a351b9"   # n     = 4  (ΟΧΙ το n_ox/n_er!)
B_ON  = 4.0                                        # b, σταθερό στο Kd sweep

STRESS = {
    "ox": np.array([10., 20., 40., 80., 150., 271., 400., 600.]),      # μM H2O2
    "er": np.array([130., 260., 500., 900., 1500., 2200., 2800., 3300.]),  # μM DTT
}

SWEEPS = [
    # (ετικέτα, param id, πλέγμα τιμών στο SBML, μετατροπή σε "τιμή παραμέτρου")
    ("K_R", P_KR, [0.022, 0.066, 0.33, 0.44, 1.32, 2.2],      lambda x: x),
    ("Kd",  P_UW, [0.01, 0.1, 0.75, 1.0, 10.0, 100.0],        lambda x: x / B_ON),
    ("n",   P_N,  [2.0, 2.5, 3.0, 3.5, 4.0],                  lambda x: x),
]

T_END, N_POINTS = 400, 4000


def hill4(S, base, amp, K, n):
    return base + amp * S**n / (K**n + S**n)


def metrics(S, y):
    """Επιστρέφει (n_eff, R2_γραμμικού, fold_change).
    ΠΡΟΣΟΧΗ: το R2 είναι του ΓΡΑΜΜΙΚΟΥ fit (μετρικό γραμμικοποίησης),
    ΟΧΙ του Hill fit — αυτό βγαίνει πάντα ~1 και δεν λέει τίποτα."""
    p0 = [y.min(), max(y.max() - y.min(), 1e-9), float(np.median(S)), 1.5]
    try:
        popt, _ = curve_fit(hill4, S, y, p0=p0, maxfev=400000)
        n_eff = abs(popt[3])
    except Exception:
        n_eff = float("nan")
    A = np.vstack([S, np.ones_like(S)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ coef
    r2 = 1.0 - resid.var() / y.var() if y.var() > 0 else float("nan")
    fc = y.max() / y.min() if y.min() > 0 else float("nan")
    return n_eff, r2, fc


def dose_response(path, sid, out_id, param_id, param_value, stress_grid):
    """Μία ολόκληρη καμπύλη. Το μοντέλο ξαναφορτώνεται σε κάθε σημείο:
    το r.reset() ΔΕΝ επαναφέρει παραμέτρους, οπότε χωρίς αυτό οι τιμές του
    προηγούμενου βήματος θα διέρρεαν σιωπηλά στο επόμενο."""
    y = []
    for S in stress_grid:
        r = roadrunner.RoadRunner(path)
        r.resetToOrigin()
        r[param_id] = param_value
        r[sid] = float(S)
        res = r.simulate(0, T_END, N_POINTS, ["time", out_id])
        y.append(res[-1, 1])
    return np.array(y)


def selftest():
    """Baseline sanity check πριν τις ~240 προσομοιώσεις."""
    print("=== SELF-TEST (baseline) ===")
    exp = {"ox": (1.686, 0.9555, 3.00), "er": (1.985, 0.9919, 2.76)}
    ok = True
    for v in ("ox", "er"):
        y = dose_response(MERGED[v], IDS[v]["S"], IDS[v]["out"], P_KR, 0.44, STRESS[v])
        ne, r2, fc = metrics(STRESS[v], y)
        e = exp[v]
        good = (abs(ne - e[0]) < 0.05 and abs(r2 - e[1]) < 0.02 and abs(fc - e[2]) < 0.10)
        ok = ok and good
        print(f"  {v}: n_eff={ne:.4f} (αναμ. {e[0]})  R2={r2:.4f} (αναμ. {e[1]})  "
              f"fold={fc:.4f} (αναμ. {e[2]})  {'OK' if good else '<<< ΑΠΟΚΛΙΣΗ'}")
        if y.max() - y.min() < 1e-6:
            print("  [!] Η έξοδος ΔΕΝ μεταβάλλεται — έλεγξε ότι οι κανόνες πέρασαν στο merge.")
            ok = False
    print(f"  → {'PASS' if ok else 'FAIL'}\n")
    return ok


def main():
    if not selftest():
        print("Το self-test απέτυχε. Δεν συνεχίζω.")
        return

    rows = []
    for variant in ("ox", "er"):
        path, sid, out_id = MERGED[variant], IDS[variant]["S"], IDS[variant]["out"]
        S = STRESS[variant]
        for label, pid, grid, to_param in SWEEPS:
            print(f"--- {variant} · {label} ---")
            for sbml_val in grid:
                y = dose_response(path, sid, out_id, pid, sbml_val, S)
                ne, r2, fc = metrics(S, y)
                rows.append({"variant": variant, "sweep": label,
                             "sbml_value": sbml_val, "param_value": to_param(sbml_val),
                             "n_eff": ne, "R2_linear": r2, "fold_change": fc,
                             "y_min": y.min(), "y_max": y.max()})
                print(f"  SBML={sbml_val:<8g} param={to_param(sbml_val):<8g}  "
                      f"n_eff={ne:7.4f}  R2={r2:7.4f}  fold={fc:7.4f}")
            print()

    with open("sweep_results.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"Γράφτηκαν {len(rows)} γραμμές στο sweep_results.csv")
    print("\nΑντιγραφή στο workbook: στήλες n_eff / R2_linear / fold_change "
          "στα κίτρινα κελιά του αντίστοιχου sweep sheet και variant block.")


if __name__ == "__main__":
    main()

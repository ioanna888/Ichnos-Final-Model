"""
ICHNOS — shared configuration (paths, constants). No dependencies on any
other ichnos module; everything else imports FROM here.

The original single-file module docstring is preserved below for context.

----------------------------------------------------------------------

ICHNOS — τρέχει τα ξεχωριστά SBML submodels ΜΑΖΙ, προγραμματιστικά.

Κάθε φορά που τρέχει:
  1. Διαβάζει τα 3 αρχεία (Ioanna, sensing module, reporter module) από το δίσκο
  2. Τα ενώνει ΣΤΗ ΜΝΗΜΗ (name-based matching, όχι id-based — βλ. GUID ιστορία)
  3. Φορτώνει το αποτέλεσμα κατευθείαν σε roadrunner και τρέχει simulation
  4. (ΝΕΟ, 2026-08-12) Σώζει επίσης το merged SBML — ΟΧΙ μόνο σαν φευγαλέο
     debug snapshot πια, γιατί το integration θα τροφοδοτήσει sensitivity
     analysis και το no-feedback comparison circuit, που χρειάζονται σταθερή
     αναφορά σε συγκεκριμένο αρχείο:
       - exportsbml\\merged_<variant>.sbml                → 'latest', ξαναγράφεται κάθε run
       - exportsbml\\archive\\merged_<variant>_<ts>.sbml    → ποτέ δεν σβήνεται/ξαναγράφεται
       - exportsbml\\archive\\merged_<variant>_<ts>.manifest.json
         → ποια source .sbml αρχεία (με mtimes) και ποιες τιμές βασικών
           παραμέτρων (mu, n, k_deg_TIP, ...) μπήκαν σε ΑΥΤΟ το merge
     Το merge πάντα ξαναχτίζεται από τα 3 πρωτότυπα αρχεία σε κάθε run· το
     script ποτέ δεν διαβάζει τίποτα πίσω από το exportsbml\\ folder. Αν δεν το θες:
     `python run_ichnos.py --no-save`.

Αν κάποιος συνάδελφος αλλάξει το δικό του .sbml αρχείο, το επόμενο run το
πιάνει αυτόματα — δεν χρειάζεται να ξαναφτιαχτεί/ξανασωθεί τίποτα ενδιάμεσο.

ΑΡΧΙΤΕΚΤΟΝΙΚΗ (CONFIRMED 2026-07-31, wet lab): Cell A (oxidative circuit:
pSynOS4 → TIP_OX → TetR → reporter) και Cell B (ER circuit: pKAR2/pUPR →
TIP_ER → TetR → reporter) είναι δύο ΞΕΧΩΡΙΣΤΑ κύτταρα/κυκλώματα — όχι δύο
ταυτόχρονα σήματα στο ίδιο κύτταρο. Δεν υπάρχει σενάριο (σκόπιμο ή ακούσιο)
όπου ένα strain βλέπει και τα δύο stress signals μαζί. Άρα οι δύο ανεξάρτητες
παραλλαγές (er/ox) παρακάτω είναι η σωστή, πλήρης αναπαράσταση — όχι
απλοποίηση εν αναμονή επιβεβαίωσης.

ΤΡΕΧΟΥΣΕΣ ΤΙΜΕΣ ΤΟΥ TIP_TETR_MODEL (CONFIRMED 2026-08-11, screenshot Ioanna
από το COPASI/Tellurium parameters panel — model scope):
    b            = 4       1/(nanomole/liter)/hour
    u_w          = 1       1/hour
    delta_w      = 0.35    1/hour ->  0.12    1/hour   (καθαρή αποικοδόμηση· +mu στην αντίδραση)
    a_TetR       = 50      nanomole/liter/hour
    k_deg_TetR   = 0.35    1/hour -> k_deg_TetR   = 0.12    1/hour          (καθαρή αποικοδόμηση· +mu στην αντίδραση)
    K_R          = 0.44    nanomole/liter
    n            = 4       dimensionless   (TetR cooperativity, βλ. πίνακα παραμέτρων)
    k_deg_TIP    = 0.30    1/hour          (ΗΔΗ combined rate: degradation+dilution) -> αφαιρέθηκε, ζει στα sensing modules (1.0 1/hour)
    P            = 1       dimensionless   (boundary/assignment rule species — ΟΧΙ constant param)
    P_min        = 0.0056  dimensionless
    mu           = 0.20    1/hour          (dilution — τώρα ΚΑΝΟΝΙΚΟ model parameter, όχι μόνο
                                             προστιθέμενο term μέσα σε reactions) ->  0.35    1/hour          (αραίωση· διπλασιασμός 120 min)

⚠️ OPEN / PENDING ITEMS:
   - Το `mu` (0.2) και το `P` (rule-governed, τιμή 1.0) είναι πλέον σωστά
     αναγνωρισμένα ως SHARED_PARAM_NAMES — κάθε submodel που ορίζει δικό του
     αντίγραφο ενοποιείται αυτόματα με το πραγματικό του Ioanna model, με
     loud warning αν οι τιμές διαφέρουν (π.χ. reporter's local P=0.5 vs
     Ioanna's πραγματικό P=1.0 — επιβεβαιωμένο conflict, σωστά χειρισμένο).
   - CONFIRMED 2026-08-12 (πραγματικό end-to-end test με τα 4 πραγματικά
     .sbml αρχεία): και οι δύο variants (er, ox) κάνουν merge, φορτώνουν σε
     roadrunner, ΚΑΙ τρέχουν simulation επιτυχώς, μηδέν dangling symbol
     references.
   - RESOLVED 2026-08-17 (Option B + cross-variant check): το `k_deg_TIP`
     είναι ο ρυθμός αποικοδόμησης του TIP ΜΕΤΑ την παραγωγή του. Η παραγωγή
     TIP γίνεται ΜΟΝΟ στα sensing modules (oxidative_module / ERModule /
     copper), ΠΟΤΕ στο binding module (TIP_TETR_MODEL) — αυτό είναι καθαρά
     δομικό/θερμοδυναμικό (binding kinetics). Άρα το `k_deg_TIP` ανήκει
     εννοιολογικά στα sensing modules και ζει ΜΟΝΟ εκεί.
       * ΑΠΟΦΑΣΗ: αφαιρείται τελείως από το TIP_TetR_binding.sbml (ήταν
         0-αναφορών/νεκρό εκεί — καμία αντίδραση αποδόμησης TIP σε εκείνο το
         αρχείο). ΔΕΝ μπαίνει στο SHARED_PARAM_NAMES (η εναλλακτική Option A
         εξετάστηκε και απορρίφθηκε: θα κρατούσε τον rate "φιλοξενούμενο" σε
         δομικό module που δεν του ανήκει).
       * ΣΥΝΕΠΕΙΑ: αφού τα sensing modules δεν κάνουν ΠΟΤΕ merge μαζί (ox/er/
         cu = ξεχωριστά κύτταρα, ξεχωριστά builds), ο ενδο-merge μηχανισμός
         [!] CONFLICT ΔΕΝ μπορεί πια να συγκρίνει το k_deg_TIP ανάμεσα σε
         variants. Η τιμή ΟΜΩΣ πρέπει να είναι ίδια παντού (ΙΔΙΟ TIP coding
         sequence σε ox/er/copper — confirmed wet lab 2026-08-17).
       * ΕΓΓΥΗΣΗ: η κοινή τιμή επιβάλλεται από το
         check_cross_variant_shared_values() — post-build validation που
         διαβάζει το k_deg_TIP από ΚΑΘΕ χτισμένο variant, τα συγκρίνει μεταξύ
         τους απευθείας, και κάνει raise αν αποκλίνουν. Αυστηρότερο από το A:
         συγκρίνει variants directly, όχι έμμεσα μέσω binding module.
       * TODO (χωριστό, δεν έγινε ακόμα): πρέπει να αφαιρεθεί χειροκίνητα το
         k_deg_TIP από το ίδιο το TIP_TetR_binding.sbml (COPASI side, Ioanna)
         — μέχρι τότε παραμένει νεκρό εκεί, αβλαβές αλλά όχι καθαρό.

Usage:
    python run_ichnos.py            # τρέχει και τις δύο παραλλαγές (er, ox), σώζει exportsbml\\*.sbml
    python run_ichnos.py er         # μόνο ER-stress variant
    python run_ichnos.py ox         # μόνο oxidative-stress variant
    python run_ichnos.py --no-save  # χωρίς να γράψει τίποτα σε exportsbml\\
"""

import os

TIP_TETR_MODEL = r"C:\Users\ioann\OneDrive - University of Patras\igem\SynBiology\integration\TIP_TetR_binding.sbml"
REPORTER_MODEL = r"C:\Users\ioann\OneDrive - University of Patras\igem\SynBiology\integration\reporter_module_v2.sbml"

TIP_TETR_MODEL_ALT_CHECK = r"C:\Users\ioann\OneDrive - University of Patras\igem\SynBiology\matlab\tiptetrbinding\TIP_TetR_binding.sbml"

VARIANTS = {
    "er": {"sensing_file": r"C:\Users\ioann\OneDrive - University of Patras\igem\SynBiology\integration\ERmodule.sbml", "tip_name": "TIP_er"},
    "ox": {"sensing_file": r"C:\Users\ioann\OneDrive - University of Patras\igem\SynBiology\integration\oxidative_module_v3.sbml", "tip_name": "TIP_ox"},
}

SHARED_PARAM_NAMES = {"mu", "P"}
SHARED_PARAM_TOLERANCE = 1e-9

EXPORT_SBML_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "exportsbml")

CROSS_VARIANT_SHARED_NAMES = {
    "k_deg_TIP",   # post-production degradation of TIP; same TIP coding sequence
                   # in ox/er/copper (confirmed wet lab 2026-08-17) → must match.
}
_CROSS_VARIANT_ABS_TOL = 1e-9

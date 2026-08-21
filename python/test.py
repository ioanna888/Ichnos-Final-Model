import tellurium as te

r = te.loadSBMLModel(r"C:\Users\ioann\OneDrive - University of Patras\igem\SynBiology\matlab\tiptetrbinding\TIP_TetR_binding_CLEAN.sbml")

print("Species IDs:", r.getFloatingSpeciesIds())
print("Parameter IDs:", r.getGlobalParameterIds())

r.reset()
result = r.simulate(0, 200, 500)
print("\nFinal values:")
for name, val in zip(result.colnames, result[-1]):
    print(f"  {name}: {val:.4f}")

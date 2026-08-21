import libsbml
import re

INPUT = r"C:\Users\ioann\OneDrive - University of Patras\igem\SynBiology\matlab\tiptetrbinding\TIP_TetR_binding.sbml"
OUTPUT = r"C:\Users\ioann\OneDrive - University of Patras\igem\SynBiology\matlab\tiptetrbinding\TIP_TetR_binding_CLEAN.sbml"

def sanitize(name):
    """Κάνει το name έγκυρο SBML SId: γράμματα/αριθμοί/_, δεν ξεκινά με αριθμό."""
    s = re.sub(r'[^A-Za-z0-9_]', '_', name)
    if re.match(r'^[0-9]', s):
        s = '_' + s
    return s

def rename_in_ast(node, mapping):
    if node is None:
        return
    if node.isName() and node.getName() in mapping:
        node.setName(mapping[node.getName()])
    for i in range(node.getNumChildren()):
        rename_in_ast(node.getChild(i), mapping)

doc = libsbml.readSBMLFromFile(INPUT)
model = doc.getModel()

id_map = {}  # old GUID id -> new clean id

for c in model.getListOfCompartments():
    new_id = sanitize(c.getName()) if c.getName() else c.getId()
    id_map[c.getId()] = new_id
for s in model.getListOfSpecies():
    id_map[s.getId()] = sanitize(s.getName())
for p in model.getListOfParameters():
    id_map[p.getId()] = sanitize(p.getName())
for r in model.getListOfReactions():
    # τα reactions δεν έχουν σημαντικό "Name" σε πολλά exports· κράτα το id αν δεν υπάρχει καλό name
    id_map[r.getId()] = r.getId()

# Ξαναγράψε compartments/species/parameters με τα νέα ids
for c in model.getListOfCompartments():
    c.setId(id_map[c.getId()])
for s in model.getListOfSpecies():
    old = s.getId()
    s.setId(id_map[old])
    s.setCompartment(id_map.get(s.getCompartment(), s.getCompartment()))
for p in model.getListOfParameters():
    p.setId(id_map[p.getId()])

# Ξαναγράψε reactions: reactants/products/kinetic law formulas
for r in model.getListOfReactions():
    for i in range(r.getNumReactants()):
        ref = r.getReactant(i)
        ref.setSpecies(id_map.get(ref.getSpecies(), ref.getSpecies()))
    for i in range(r.getNumProducts()):
        ref = r.getProduct(i)
        ref.setSpecies(id_map.get(ref.getSpecies(), ref.getSpecies()))
    kl = r.getKineticLaw()
    if kl is not None:
        rename_in_ast(kl.getMath(), id_map)

# Ξαναγράψε Rules
for rule in model.getListOfRules():
    var = rule.getVariable()
    if var in id_map:
        rule.setVariable(id_map[var])
    rename_in_ast(rule.getMath(), id_map)

libsbml.writeSBMLToFile(doc, OUTPUT)
print("File written:", OUTPUT)
print("\nChecking — new species ids:", [s.getId() for s in model.getListOfSpecies()])
print("New Parameter ids:", [p.getId() for p in model.getListOfParameters()])


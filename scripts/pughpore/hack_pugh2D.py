# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 11:42:54 2016

@author: gregor
"""
import nanopores.tools.box as box
import nanopores.tools.balls as balls
import nanopores.geometries.pughpore as pugh
import nanopores.py4gmsh as gmsh
Box = balls.Box

pugh.set_tol(1e-5)

dom = pugh.get_domain_cyl()
dom.addboundaries(leftb=dom.boundary("left"))
left = dom.getboundary("leftb")

mol = pugh.EmptySet()
dom.addsubdomain(mol, "molecule")
dom.addboundary(mol.boundary() - left, "moleculeb")

dom.compute_entities()
dom.compute_boundaries(True)

# points of half circle
x0 = [0.,-15.]
r = 2.
lcMolecule = 0.4
lc = 1.

def entity2box(ent):
    intervals = [(f if isinstance(f, tuple) else (f,f)) for f in ent]
    return Box(intervals=intervals)
    
edgeinds = list(left.indexsets[1])
edgeents = [dom.entities[1][i] for i in edgeinds]
edge = [entity2box(dom.entities[1][i]) for i in edgeinds]
points = [(x0[0], x0[1]-r), tuple(x0), (x0[0], x0[1]+r)]
circle = [Box(points[i], points[i+1]) for i in range(len(points)-1)]
N = len(edge)

dic = box.multi_box_union(edge + circle)

# add additional point entities
for p in dic["entities"][0]:
    if not p in dom.entities[0]:
        dom.entities[0].append(p)

# add new edge entities and compute replacement
replace = {i:[] for i in edgeinds}
circleb = []
for s, ent in zip(dic["esets"][1], dic["entities"][1]):
    for j in s:
        if j < len(edgeinds): # is old edge
            i = edgeinds[j]
            replace[i].append(ent)
        if j >= len(edgeinds): # belongs to circle
            print j
            circleb.append(ent)
for k in replace.keys():
    for i, ent in enumerate(replace[k]):
        if ent in dom.entities[1]:
            j = dom.entities[1].index(ent)
        else:
            dom.entities[1].append(ent)
            j = len(dom.entities[1]) - 1
            print j, ent
        replace[k][i] = j
for k, v in replace.items():
    if len(v)==1 and k==v[0]:
        replace.pop(k)
print replace
old = set(replace.keys())
new = box.union(set(v) for v in replace.values())
# replace edge indices in boundary
left.indexsets[1] = left.indexsets[1] - old | new

# compute left circle boundary
for i, ent in enumerate(circleb):
    if ent in dom.entities[1]:
        j = dom.entities[1].index(ent)
    else:
        dom.entities[1].append(ent)
        j = len(dom.entities[1]) - 1
    circleb[i] = j
print "circle:", circleb

# gmsh circle
lcCirc = lcMolecule*lc
m0, m1 = x0[0], x0[1]
pcirc = [(m0, m1), (m0, m1-r), (m0+r, m1), (m0, m1+r)]
dom.entities[0].append(pcirc[2])

dom.gmsh_entities = [[None for e in k] for k in dom.entities]
pcirc = [dom.entity_to_gmsh(p, 0, lcCirc) for p in pcirc]

surfs = [gmsh.Circle([pcirc[1], pcirc[0], pcirc[2]]),
         gmsh.Circle([pcirc[2], pcirc[0], pcirc[3]])]
dom.gmsh_entities[1] += surfs
N = len(dom.gmsh_entities[1])
circlearc = [N-2, N-1]

for k, v in replace.items():
    removed = False
    for j in list(v):
        print "adding", j,"to replace"
        if j in circleb:
            replace[k].remove(j)
            removed = True
    if removed:
        replace[k].extend(circlearc)
        
print replace
# add edge indices to molecule boundary
mol.bdry().indexset = set(circleb + circlearc)
mol.bdry().indexsets[1] = set(circleb + circlearc)
for i in circleb:
    mol.bdry().orients[i] = -1
for i in circlearc:
    mol.bdry().orients[i] = 1

# replace edge indices sub.boundaries
for sub in dom.subdomains + dom.boundarysubs:
    iset = sub.bdry().indexset
    orients = sub.bdry().orients
    for i in iset & old:
        iset.remove(i)
        for j in replace[i]:
            iset.add(j)
            if j in circlearc:
                orients[j] = -1
            else:
                orients[j] = orients[i]
            print sub.name, i, j, orients[j]

dom.entities_to_gmsh_merge(lc)
# rebuild boundaries involving balls
for bou in dom.boundaries:
    bou.indexset = bou.csg.evalsets()[1]
    
dom.physical_to_gmsh(True)  
dom.geo = box.to_mesh()
dom.geo.params = dom.params
if hasattr(dom, "synonymes"):
    dom.geo.import_synonymes(dom.synonymes) 
dom.plot()






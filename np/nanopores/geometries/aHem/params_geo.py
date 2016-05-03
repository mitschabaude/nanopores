""" --- geometry parameters for Howorka cylindrical geometry --- """

nm = 1e0

# @TODO maybe include tolc in parent file
tolc = 1e-14*nm  # tolerance for coordinate comparisons

dim = 3

synonymes = {
    "pore": {"poretop", "porecenter", "porebottom"},
    "bulkfluid": {"fluid_bulk_top", "fluid_bulk_bottom"},
    "fluid": {"pore", "bulkfluid"},
    "solid": {"membrane", "ahem", "molecule"},
    "protein": "ahem",
    "proteinb": "ahemb",
    "noslip": {"ahemb", "membraneb", "moleculeb"},
    "bulk": {"upperb", "lowerb"},
    "nopressure": "bulk",
    "ground": "upperb",
    "bV": "lowerb",
    "ions": "fluid",
    "lipid": "membrane",
    "exittime": "fluid",
    "exit": "poreexit",
    "sideb": {"uppersideb", "lowersideb"},
    "upperbulkb": {"upperb", "uppersideb"},
    "lowerbulkb": {"lowerb", "lowersideb"},
}

# characteristic length / mesh size h
lc = nm
lcMolecule = lc*3e-2
lcCenter = lc*5e-2
lcOuter = lc*1e-1
lcMembrane = lc*5e-3

# provide default values for boundary layer around membrane/molecule
membraneblayer = None
moleculeblayer = None

#effective pore radius

'''
Sketch of alpha-Hemolysin domain with dimensions
____________________________________________ 
|                                  .        |           
|                                  .        |           
/---------------------R------------.--------/           
|                                  .        |           
|                                  .        |           
|                                 l3        |           
|                                  .        |           
|                                  .        |           
/--------------r3--------------/   .        |           
|             _________________ ....        |
|            |    .            |   .        |           
|            |    .            |   .        |           
|            |    .            |   .        |           
|            |   l0            |   .        |           
/-----r1-----/    .            |  l2        |             
|            |    .            |   .        |             
|            |    .            |   .        |             
|       _____|.....            |   .        |                  
|       |         .            |   .        |             
|       |         .       _____|....        |           
|       |         .      |                  |                  
/---r0--/         .      |__________________|...                            
|       |        l1                            lm       
|       |         .       __________________ ...                 
|       |         .      |                  |           
|       |_________.______|..........        |           
|                                  .        |           
/-----------r2-----------/        l4        |           
|__________________________________.________|

EDGES:

  .-f--.
  .    .___________
  .    /
  .   /
  .  /
  . /
  ./
  |
  |
'''
f = 0.5*nm

l0 = 4.0*nm

l1 = 5.5*nm

l2 = 6.2*nm

lm = 2.2*nm

l3 = 15.0*nm

l4 = 1.0*nm

r0 = 1.1*nm

r1 = 1.3*nm

r2 = 1.7*nm

r3 = 4.2*nm

R = 15.0*nm

rMolecule = 2.2*nm

# (c) 2016 Gregor Mitscha-Baude
"back-of-the-envelope computation of current through pughpore"
import nanopores
phys = nanopores.Physics("electrolyte")

# parameters
L = 9. #35.
diam = 2. #7.5

nano = 1e-9
L = nano*L # m
bV = 0.02 # V
#A = nano**2 * diam**2 # m**2
A = nano**2 * diam**2 * 3.141592/4 # m**2
c0 = 1000 # mol/m**3
mol = phys.mol # 1/mol

C0 = c0*mol # 1/m**3
q = phys.qq # C
kT = phys.kT # J = m**2 / s**2 * kg = N*m
D = 0.5*(1.96 + 2.03) * 1e-9 # m**2/s
E = 0.9* bV / L # V/m

mu = D / kT # (m**2 / s) / J = (m/s) / N
v = mu * q*E # (m/s) / N * N = m/s

j = 2.*C0*v # 1/m**3 * m/s = 1/(m**2 s)
J = q*j*A # C/s = A

print "current:", J*1e12, "pA"
print "conductance:", J/bV*1e9, "nS"
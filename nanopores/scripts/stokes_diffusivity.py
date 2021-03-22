''' theoretical diffusivity from stokes law compared with measurements '''

from math import pi

k = 1.38064e-23
T = 273 + 20
eta = 1e-3

def Dstokes(r):
    return k*T/(6*pi*eta*r) * 1e18
def rstokes(D):
    return k*T/(6*pi*eta*D) * 1e18

# experimental data
# name : radius [nm], diffusivity [nm**2/ns]
data = {
    "K+" : (0.152, 1.96),
    "Na+" : (0.116, 1.33),
    "Cl-" : (0.167, 2.03),
    "ion" : (0.11, 1.9)
}

if __name__ == "__main__":
    # some examples

    s0 = "%s\nD (measured): %s\nD (Stokes): %s"
    s1 = "r (measured): %s\nr (Stokes): %s\n"
    for name, (r, D) in list(data.items()):
        Ds = Dstokes(r)
        print(s0 % (name, D, Ds))
        rs = rstokes(D)
        print(s1 % (r, rs))

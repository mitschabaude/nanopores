""" --- geometry parameters for nanowire geometry --- """

nm = 1e-9

# tolerance for coordinate comparisons
tolc = 1e-10*nm

# box width (x coordinate)
w = 176*nm
# nanowire length (y coordinate)
l = 1000*nm
# box height (z coordinate)
h = 253*nm

# oxide height below nanowire
h_side = 145*nm
# oxide width left and right of nanowire
w_side = 50*nm

# dimensions of oxide layer on nanowire
w_layer = 8*nm
h_layer = 8*nm

# dimensions of silicon core of nanowire
w_core = 60*nm
h_core = 50*nm

# position of silicon core
x_core = w_side + w_layer
z_core = h_side

# total dim of nanowire
w_wire = w_core + 2*w_layer
h_wire = h_core + h_layer

# TODO: approach reasonable?
# effective dopant radius
# chosen so that volume of one regular tetraeder of this side length equals
# volume of ball of actual radius
pi = 3.141592653589793
#r_dop = 0.5*nm
#r_eff = r_dop*2*(2**.5*pi)**(1./3) #TODO
r_eff = 3*nm

# characteristic mesh width scales
lc = 16*nm
rellcwire = 1./8
#lcwire = r_eff


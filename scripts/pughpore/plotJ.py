# (c) 2016 Gregor Mitscha-Baude
import nanopores.tools.fields as f
import matplotlib.pyplot as plot

field = f.get_fields("pughcenter")
z = [x[2] for x in field["x"]]
J = [j*1e12 for j in field["J"]]
I = sorted(range(len(z)), key=lambda k: z[k])
z1 = [z[i] for i in I]
J1 = [J[i] for i in I]

plot.plot(z1, J1, "s-")
plot.xlabel("z position of molecule [nm]")
plot.ylabel("current [pA]")
plot.title("current at -100mV for molecule along pore center")
plot.show()

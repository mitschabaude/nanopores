# (c) 2016 Gregor Mitscha-Baude
import nanopores.tools.fields as f
import matplotlib.pyplot as plot

field = f.get_fields("pughcenter", bulkcon=1e3)
z = [x[2] for x in field["x"]]
J = [j*1e12 for j in field["J"]]
I = sorted(range(len(z)), key=lambda k: z[k])
z1 = [z[i] for i in I]
J1 = [J[i] for i in I]

plot.plot(z1, J1, "s-")
plot.xlabel("z position of molecule [nm]")
plot.ylabel("current [pA]")
plot.title("current at -100mV for molecule along pore center")

for h in [2., 3., 4.]:
    Jopen = f.get_fields("pughopen", bulkcon=1e3, h=h)["J"][0]
    plot.plot(z1, [Jopen*1e12]*len(z1), "--", label="open, h=%s" %h)
plot.legend()
plot.show()

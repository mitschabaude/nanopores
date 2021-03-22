# (c) 2016 Gregor Mitscha-Baude
import nanopores.tools.fields as f
import matplotlib.pyplot as plot
import folders

f.update()
for h in [2.,4.]:
    field = f.get_fields("pughcenter", bulkcon=1e3, Qmol=8, h=h)
    z = [x[2] for x in field["x"]]
    J = [j*1e12 for j in field["J"]]
    I = sorted(list(range(len(z))), key=lambda k: z[k])
    z1 = [z[i] for i in I]
    J1 = [J[i] for i in I]
    for i in I:
        print(z[i], J[i])
    
    plot.plot(z1, J1, "s-", label="molecule, h=%s" %h)
    plot.xlabel("z position of molecule [nm]")
    plot.ylabel("current [pA]")
    plot.title("current at -100mV for molecule along pore center")

for h in []: #[1., 2., 3., 4.]:
    Jopen = f.get_fields("pughopen", Qmol=8, h=h)["J"][0]
    plot.plot(z1, [2.*Jopen*1e12]*len(z1), "--", label="open, h=%s" %h)
    
plot.legend()
plot.show()

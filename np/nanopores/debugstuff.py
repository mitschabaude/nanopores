from dolfin import *

__all__ = ["debugform", "printfunc", "searchinproblem"]

def debugform(a):
	#print str(a.function_spaces()) hat nur dolfin.Form
	d = a.subdomain_data()
	#print str(d.values()[0])
	#for itg in a.integrals():
	#	print "Integral: %r, %r, %r, %r, %r" % (itg._integral_type,
        #itg._domain, itg._subdomain_id, itg._metadata, itg._subdomain_data)
	print len(d[list(d)[0]])
	for argument in a.arguments():
		print str(argument),", ",type(argument),", ", \
              str(argument.number()),", ", str(argument.part()), ", ", \
              str(argument.function_space()) 
	for coeff in a.coefficients():
		if(isinstance(coeff, Function)) : 
			print "This is a Function:"
			print str(coeff.name()),", ",type(coeff),", ",\
                str(coeff.function_space())
		else: print str(coeff.name()),", ",type(coeff),", ", str(coeff.element())
	for key,value in a.subdomain_data().items():
		print repr(key),": "
		for key2,value2 in value.items():
			#plot(value2,title=repr(key2),interactive=True)
			print str(key2),": ",str(value2.id())#.mesh().id())
	#assemble(a)

    
def printfunc(f):
	print f.name(),": dim",f.function_space().dim(),"id",f.id(),"count",f.count()
	#for i in range(f.function_space().num_sub_spaces()):
	#	printfunc(f.sub(i))

def searchinproblem(f,problem):
	for c in set(problem.a.coefficients() + problem.L.coefficients()):
		print c.name(),": count", c.count()
		if (c.count() == f.count()):
			print "matching count found in bilinear form"


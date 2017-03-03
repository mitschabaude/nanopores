import matplotlib
matplotlib.use("Agg")
import nanopores as nano
import nanopores.geometries.pughpore as pughpore
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import nanopores.tools.fields as f
HOME = os.path.expanduser("~")
PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
FIGDIR = os.path.join(PAPERDIR, "figures", "")
DATADIR = os.path.join(HOME,"Dropbox", "nanopores", "fields")
f.set_dir(DATADIR)

number=False

geop = nano.Params(pughpore.params)
hpore=geop.hpore
#params=dict(avgbind=8.7e6,P_bind=5.e-3,z0=hpore/2.+5.)
def save_fig(params,fieldsname):
	figname = fieldsname+'_%.1e_%.1e_%.1e_%.1e_'%(params["avgbind1"],params["avgbind2"],params["P_bind1"],params["P_bind2"],)+str(params["z0"])+'.eps'
	data=f.get_fields(fieldsname,**params)
	t1 = data["t1"]
	a1 = data["a1"]
	t2 = data["t2"]
	a2 = data["a2"]
	print len(t1)
	print len(t2)
#	print 't1'
#	print t1
#	print 'a1'
#	print a1
#	print 't2'
#	print t2
#	print 'a2'
#	print a2
#	t1 = np.array([])
#	t2 = np.array([])
#	a1 = np.array([])
#	a2 = np.array([])
#	print '#data = %i'%len(T)
#	for i in range(len(T)):
#	    T_=np.array(T[i])*1e-6
#	    J_=np.array(J[i])
#	    tau_off=np.sum(T_)
#	    amp = (2060.-np.inner(J_,T_)/tau_off)/2060.*100
#	    if number: ax.text(tau_off,amp,'%i'%i,fontsize=9)
#	    if tau_off<.1:
#		t1=np.append(t1,np.array([tau_off]))
#		a1=np.append(a1,np.array([amp]))
#	    else:
#		t2=np.append(t2,np.array([tau_off]))
#		a2=np.append(a2,np.array([amp]))



	plt.plot([5e-4,.7],[0.,0.],linewidth=2,color='lightgreen')
	plt.plot([1.,1e2],[0.,0.],linewidth=2,color='green')
	ax=plt.gca()
	plt.scatter(t1,a1,color='lightgreen')
	plt.scatter(t2,a2,color='green')
#	ax.set_xlim([1e-4,500.])
#	ax.set_ylim([-0.40,4.0])
	ax.set_xscale('log')
	ax.invert_yaxis()
	ax.set_xlabel(r'$\tau_{off}$ [ms]',fontsize=15)
	ax.set_ylabel(r'A/I$_0$ [%]',fontsize=15)
	ax.text(.011,-0.03,'I',fontsize=15)
	ax.text(5.,-0.03,'II',fontsize=15)
	plt.tight_layout()
	nano.savefigs(name=figname,DIR='/home/lv70496/benjamin/plots/')
	print 'savefig:'
	print figname

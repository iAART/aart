import numpy as np
import subprocess

pythoncmd="python "
cmd0="ChangeParams.py --p_image=0 --limits=100 --dx0=0.02 --dx1=20 --dx2=20"
subprocess.run([pythoncmd+cmd0],shell=True)

cmd2="lensingbands.py"
cmd3="raytracing.py"
cmd4="gfactor.py"
cmd5="Radius_Redshifit.py"
cmd6="rm -r Results/*"

spinvals=[0.02,0.2,0.5,0.7,0.8,0.85,0.9,0.94,0.98,0.99]
incvals=[20,30,40,50,60,70,75,80,85,89]

#spinvals=[0.02,0.2,0.4,0.6,0.8,0.9,0.95,0.99]
#incvals=[40,50,60,70,80,85]

results=[]

for incval in incvals:
	for spinval in spinvals:

		if incval>40:
			cmd0="ChangeParams.py --p_image=0 --limits=%s --dx0=%s --dx1=20 --dx2=20"%(50,0.01)
			subprocess.run([pythoncmd+cmd0],shell=True)
		elif incval>60:
			cmd0="ChangeParams.py --p_image=0 --limits=%s --dx0=%s --dx1=20 --dx2=20"%(30,0.005)
			subprocess.run([pythoncmd+cmd0],shell=True)
		elif incval>70:
			cmd0="ChangeParams.py --p_image=0 --limits=%s --dx0=%s --dx1=20 --dx2=20"%(8,0.002)
			subprocess.run([pythoncmd+cmd0],shell=True)

		cmd1="ChangeParams.py --a=%s --i=%s"%(spinval,incval)

		subprocess.run([pythoncmd+cmd1],shell=True)

		subprocess.run([pythoncmd+cmd2],shell=True)

		subprocess.run([pythoncmd+cmd3],shell=True)

		subprocess.run([pythoncmd+cmd4],shell=True)

		subprocess.run([pythoncmd+cmd5],shell=True)

		results.append([spinval,incval,np.load("gmin.npy"),np.load("gmax.npy"),np.load("rmax.npy"),np.load("grmax.npy")])

		subprocess.run([cmd6],shell=True)

np.save("aig.npy",np.array(results))

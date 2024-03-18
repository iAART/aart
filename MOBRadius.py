import numpy as np
import subprocess

pythoncmd="python "
cmd0="ChangeParams.py --p_image=0 --limits=10 --dx0=0.02 --dx1=5 --dx2=5"
subprocess.run([pythoncmd+cmd0],shell=True)

cmd2="lensingbands.py"
cmd3="raytracing.py"
cmd4="gfactor.py"
cmd5="Radius_Redshifit.py"
cmd6="rm -r Results/*"

spinvals=[0.05,0.5,0.9]
incvals=[1,20,30,40,50,60,70]
results=[]

for spinval in spinvals:
	for incval in incvals:

		cmd1="ChangeParams.py --a=%s --i=%s"%(spinval,incval)

		subprocess.run([pythoncmd+cmd1],shell=True)

		subprocess.run([pythoncmd+cmd2],shell=True)

		subprocess.run([pythoncmd+cmd3],shell=True)

		subprocess.run([pythoncmd+cmd4],shell=True)

		subprocess.run([pythoncmd+cmd5],shell=True)

		results.append([spinval,incval,np.load("result.npy")])

		subprocess.run([cmd6],shell=True)

np.save("aig.npy",np.array(results))

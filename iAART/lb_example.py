import matplotlib.pyplot as plt
import numpy as np
import lb_mod as lb

lbs=lb.lb(dx0=.5,dx1=1,dx2=1,a=0.94,iobs=0.3)

plt.plot(lbs[0],lbs[1],color="k",linewidth=0.3,linestyle="--")
plt.plot(lbs[0],-lbs[1],color="k",linewidth=0.3,linestyle="--")

plt.plot(lbs[2][:,0],lbs[2][:,1],'k',linewidth=0.2)

plt.plot(lbs[3][:,0],lbs[3][:,1],'b',linewidth=0.2)
plt.plot(lbs[4][:,0],lbs[4][:,1],'b',linewidth=0.2)

plt.plot(lbs[5][:,0],lbs[5][:,1],'r',linewidth=0.2)
plt.plot(lbs[6][:,0],lbs[6][:,1],'r',linewidth=0.2)

#plt.plot(lbs[2],lbs[3],"r.")

plt.xlabel(r"$\alpha$"+" "+"(M)")
plt.ylabel(r"$\beta$"+" "+"(M)")

plt.show()
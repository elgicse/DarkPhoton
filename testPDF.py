import matplotlib.pyplot as plt

from proton_bremsstrahlung import *
from settings import * # define here the paraphoton mass and epsilon

if 'epsilon' not in vars():
	epsilon = 0.0002
if 'mDarkPhoton' not in vars():
	mDarkPhoton = 0.07

theta = math.pi / 6. # test paraphoton angle
mom = range(1,100,1) # test paraphoton momenta
dn = []

# compute differential production rate
for m in mom:
    dn.append(dNdPdTheta(m,theta,mDarkPhoton,epsilon))

# plot
plt.ion()

plt.semilogy(mom,dn)
plt.ylabel(r'$dN / dP\ d\theta$')
plt.xlabel(r'$P_A$ (GeV)')
plt.text(40, pow(10.,-10.), r'$\theta=\pi/6,\ m_A=0.07$ GeV$,\ \epsilon=0.0002$')
plt.draw()


thetas = np.linspace(-1*math.pi,1.*math.pi,360,endpoint=False).tolist()
mm = 0.01
dnt = []

for th in thetas:
	dnt.append(dNdPdTheta(mm,th,mDarkPhoton,epsilon))

plt.figure()
plt.plot(thetas,dnt)
plt.xlabel(r'$\theta_A$ (rad)')
plt.ylabel(r'$dN / dP\ d\theta$')
plt.text(-3.9, 0.2*pow(10.,-10.), r'p = 10 MeV, $m_A=0.07$ GeV$,\ \epsilon=0.0002$')
plt.draw()

#plt.show()
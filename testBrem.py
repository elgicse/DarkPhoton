from array import array
from proton_bremsstrahlung import *

mDarkPhoton = 0.5
epsilon = 1.e-6
pEn = np.linspace(100., 10000., 50)
prodRates = []
for en in pEn:
    protonEnergy = en
    protonMomentum = momentum(protonEnergy,mProton)
    pr = prodRate(mDarkPhoton, epsilon)
    print 'Proton energy is ', en, ', prod rate is ', pr
    prodRates.append(pr)

    


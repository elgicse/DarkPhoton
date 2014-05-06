import math
from proton_bremsstrahlung import *
from paraphotonDecay import *
from dark_photon import *
from settings import *
import os

#epsPower = [-3., -4., -5., -6., -7.]
#massPower = [-2., -1., 0.]
output = []
valEps =  []
valMass = []

ngen = 3500

# 0
#mEpsPairs = [[3.e-3, 1.e-3], [9.e-3, 3.e-4], [9.e-3, 2.e-7], [3.e-2, 1.e-4], [8.e-2, 3.e-5], [8.e-2, 7.e-8], [2.e-1, 1.e-5], [2.e-1, 4.e-8], [7.e-1, 1.e-6], [7.e-1, 5.e-8], [9.e-1, 1.e-7]]
# 1
#mEpsPairs = [[3.e-3, 9.e-4], [9.e-3, 2.e-4], [9.e-3, 3.e-7], [3.e-2, 9.e-5], [3.e-2, 1.e-7], [3.e-3, 3.e-7], [8.e-2, 2.e-5], [2.e-1, 9.e-6], [2.e-1, 5.e-8], [7.e-1, 9.e-7], [7.e-1, 6.e-8]]
# 2
#mEpsPairs = [[1.0, 1.e-7], [1.0, 2.e-7], [2.e-1, 1.e-5], [2.e-1, 5.e-8], [8.e-1, 7.e-8], [9.e-1, 4.e-7], [3.e-3, 8.e-4], [3.e-2, 8.e-5], [2.e-1, 8.e-6], [7.e-1, 8.e-7], [1.e-3, 3.e-7]]
# 3
#mEpsPairs = [ [3.e-3, 2.e-7], [1.e-2, 1.e-7], [4.e-2, 8.e-8], [4.e-1, 4.e-7], [6.e-1, 1.e-6], [7.e-2, 2.e-5], [1.e-1,1.5e-5] ]
# 4
#mEpsPairs = [ [2.e-3, 1.e-3], [3.e-3, 7.e-4], [6.e-3, 3.e-4], [3.e-2, 7.e-5], [1.e-1, 1.e-5], [4.e-1, 2.e-6], [2., 2.e-7] ]
# 5
#mEpsPairs = [ [4.e-1, 4.e-8], [1.e-2, 1.5e-7], [1.e-1, 7.e-6], [7.e-1, 1.5e-6], [1.e-2, 1.e-4], [3.e-2, 3.e-5], [3.e-1, 3.e-6] ]
# For the Tracking
mEpsPairs = [[8.e-1, 3.e-7]]

with open("sensitivities_1.225e7decays-6.txt","a") as ofile:

	ofile.write("Eps \t Mass \t ProdRate \t BR(ee) \t BR(mumu) \t Acc(ee) \t Acc(mumu) \t TotalAcc \t nEventsInSHiP\n")

	for mEps in mEpsPairs:
	#for mp in massPower:
		# Generate A' -> ll decays in the rest frame
		val_mDarkPhoton = mEps[0]
		val_epsilon = mEps[1]
		print "Making rest frame decays for M = %s ..."%val_mDarkPhoton

		if not os.path.isfile("out/NTuples/childrenRestFrame_e_m%s.root"%val_mDarkPhoton):
			makeNtupleDecayRestFrame(e,val_mDarkPhoton,ngen)
		if val_mDarkPhoton > 2.*(mmu/1000.):
			if not os.path.isfile("out/NTuples/childrenRestFrame_mu_m%s.root"%val_mDarkPhoton):
				makeNtupleDecayRestFrame(mu,val_mDarkPhoton,ngen)
	
		#for epsp in epsPower:
		#	val_epsilon = 1.*pow(10.,epsp)
		# total production rate of A'
		val_norm = prodRate(val_mDarkPhoton,val_epsilon)
		# number of A' produced
		val_numDarkPhotons = int(math.floor(val_norm*protonFlux))

		if not os.path.isfile("out/NTuples/ParaPhoton_eps%s_m%s.root"%(val_epsilon,val_mDarkPhoton)):
			print "Creating four-momenta for eps, mass = ",val_epsilon, val_mDarkPhoton
			create4Momenta(val_mDarkPhoton,val_epsilon,val_norm,ngen)

		# Boost children according to mother momentum
		print "Boosting children according to mother momentum..."
		#if not os.path.isfile("out/NTuples/childrenLabFrame_e_eps%s_m%s.root"%(val_epsilon,val_mDarkPhoton)):
		acc_e = boostChildren(e,val_mDarkPhoton,val_epsilon)
		#if not os.path.isfile("out/NTuples/childrenLabFrame_mu_eps%s_m%s.root"%(val_epsilon,val_mDarkPhoton)):
		if val_mDarkPhoton > 2.*(mmu/1000.):
			acc_mu = boostChildren(mu,val_mDarkPhoton,val_epsilon)
		else:
			acc_mu = 0.

		bre = leptonicBranchingRatio(val_mDarkPhoton, val_epsilon, e)
		brmu = leptonicBranchingRatio(val_mDarkPhoton, val_epsilon, mu)

		tot_acc = acc_e*bre + acc_mu*brmu

		sens = protonFlux * val_norm * tot_acc

		ofile.write("%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n"%(val_epsilon,val_mDarkPhoton,
				val_norm,bre,brmu,acc_e,acc_mu,tot_acc,sens))

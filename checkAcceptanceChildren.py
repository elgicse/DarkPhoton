# experiment settings
firstVolume = [60., 100., 2.5] # start, end, radius (m)
secondVolume = [110., 140., 2.5] # start, end, radius (m)


def inAcceptance(vtx, pChild1, pChild2):
	# Input: vertex (TVector3), children 4-momenta (TLorenzVector)
	# Check if A' -> l+l- vertex is in the first volume
	if (vtx.Z() > firstVolume[0]) and (vtx.Z() < firstVolume[1]):
		if (vtx.X()**2. + vtx.Y()**2.) < firstVolume[2]:
			# Check if child 1 goes through the detector:
			tx1 = pChild1.Px() / pChild1.Pz()
			ty1 = pChild1.Py() / pChild1.Pz()
			endPos1 = r.TVector3()
			endPos1.SetZ(firstVolume[1])
			endPos1.SetX( vtx.X() + tx1*(endPos1.Z() - vtx.Z()) )
			endPos1.SetY( vtx.Y() + ty1*(endPos1.Z() - vtx.Z()) )
			if (endPos1.X()**2. + endPos1.Y()**2.) < firstVolume[2]:
				# Check if child 2 goes through the detector:
				tx2 = pChild1.Px() / pChild1.Pz()
				ty2 = pChild1.Py() / pChild1.Pz()
				endPos2 = r.TVector3()
				endPos2.SetZ(firstVolume[1])
				endPos2.SetX( vtx.X() + tx2*(endPos2.Z() - vtx.Z()) )
				endPos2.SetY( vtx.Y() + ty2*(endPos2.Z() - vtx.Z()) )
				if (endPos2.X()**2. + endPos2.Y()**2.) < firstVolume[2]:
					return True
	# Check if A' -> l+l- vertex is in the second volume
	elif (vtx.Z() > secondVolume[0]) and (vtx.Z() < secondVolume[1]):
		if (vtx.X()**2. + vtx.Y()**2.) < secondVolume[2]:
			# Check if child 1 goes through the detector:
			tx1 = pChild1.Px() / pChild1.Pz()
			ty1 = pChild1.Py() / pChild1.Pz()
			endPos1 = r.TVector3()
			endPos1.SetZ(secondVolume[1])
			endPos1.SetX( vtx.X() + tx1*(endPos1.Z() - vtx.Z()) )
			endPos1.SetY( vtx.Y() + ty1*(endPos1.Z() - vtx.Z()) )
			if (endPos1.X()**2. + endPos1.Y()**2.) < secondVolume[2]:
				# Check if child 2 goes through the detector:
				tx2 = pChild1.Px() / pChild1.Pz()
				ty2 = pChild1.Py() / pChild1.Pz()
				endPos2 = r.TVector3()
				endPos2.SetZ(secondVolume[1])
				endPos2.SetX( vtx.X() + tx2*(endPos2.Z() - vtx.Z()) )
				endPos2.SetY( vtx.Y() + ty2*(endPos2.Z() - vtx.Z()) )
				if (endPos2.X()**2. + endPos2.Y()**2.) < secondVolume[2]:
					return True
	# Otherwise
	return False
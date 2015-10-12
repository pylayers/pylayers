import numpy as np
def SmoothMeshLine( plines, max_res, ratio=1.3,check=True,allowed_max_ratio=1.25):
	"""
	Parameters
	----------

	plines : np.array

	Returns
	-------

	lines : np.array

	"""

	dlines = np.diff(plines)
	Npoints = np.ceil(dlines/max_res)
	Npoints = Npoints.astype(int)

	for k,N in enumerate(Npoints):
		if k!=len(dlines)-1:
			l = np.linspace(plines[k],plines[k+1],N,endpoint=False)
		else:
			l = np.linspace(plines[k],plines[k+1],N+1,endpoint=True)

		try:
			lines = np.hstack((lines,l))
		except:
			lines = l
	#max_ratio = ratio*allowed_max_ratio


	if check:
		EC,pos,Etype=CheckMesh(lines,0,max_res,max_ratio,0);

	return lines

def CheckMesh(lines, min_res, max_res, ratio, verbose=False):
	""" Check if mesh lines are valid

   Parameters
   ----------

   lines   : np.array()
   min_res: minimal allowed mesh-diff
   max_res: maximal allowed mesh-diff
   ratio:   maximal allowed mesh-diff ratio
   be_quiet: disable warnings

   Returns
   -------

   EC:     error code (number of found errors)
   pos:    line positions with error
   E_type: error type

   From Matlab code of Thorsten Liebig

   See Also
   --------

	"""


	diff_lines = np.diff(lines)
	EC = 0
	E_type = []

	pos = []
	max_err = np.where(diff_lines>max_res)[0]

	if (len(max_err) >0 & verbose):
		print('CheckMesh : found resolution larger than max_res')

		pos = np.hstack((pos,max_err))
		EC = EC + len(max_err)
		E_type.append(1)

	min_err = np.where(diff_lines<min_res)[0]
	if (len(min_err)>0 & verbose):
		warning('CheckMesh : found resolution smaller than min_res')
		pos = np.hstack((pos,min_err))
		EC = EC + len(min_err)
		E_type.append(2)

	r = diff_lines[1:]/diff_lines[0:-1]

	if (r>ratio*1.01).any():
		u = np.where(r>ratio*1.01)
		pos = np.hstack((pos,u))
		E_type.append(3)

	if (r<(1/ratio)*1.01).any():
		u = np.where(r<(1/ratio)*1.01)
		pos = np.hstack((pos,u))
		E_type.append(4)


	return EC,pos,E_type

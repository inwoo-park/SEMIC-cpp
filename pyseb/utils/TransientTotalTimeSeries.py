#!/usr/bin/env python3
import numpy as np
import multiprocessing, tqdm, functools, os
from ..geometry.GetAreas import GetAreas
def AreaSum(i, data, elements, flags, Areas):
   out = np.sum(Areas*np.mean(data[elements[flags,:],i]))
   return out

def TransientTotalTimeSeries(elements, x, y, data, **kwargs):
	'''TransientTotalTimeSeries - compute total data evolution through time

	Usage
	 output = TransientTotalTimeSeries(elements,x,y,[md.results.TransientSolution.Vel]);
	 output = TransientTotalTimeSeries(elements,x,y,[md.results.TransientSolution.Vel],md.mask.ice_levelset < 0);

	# others..
	output = TransientTotalTimeSeries(md.mesh.elements,md.mesh.x,md.mesh.y,md.smb.dailydlradiation,'mask',md.mask.ice_levelset < 0,'mean',1);

	Options
	 nargin >= 6
	 mean    - mean value with weights from area. (default: 0)
	 mask    - masking for specific region. (default: md.mask.ice_levelset < 0);
	 verbose - show process for debugging (defaut: 1)
	'''

	# data should be set as npts x ntimes
	s    = np.shape(data)
	npts = len(x)
	if (npts == s[1]): # {{{
		data = np.transpose(data)
		s = np.shape(data)
	# }}}

	# check nvertices+1 for timeseries data
	if (npts+1 == s[0]): # {{{
		data = data[:-1,:]
	# }}}

	# check time length 
	if s[1] == 1: # {{{
		raise Exception('check time series of input data set. it should be more than 1');
	# }}}

	# check input argument.
	isparallel = 0
	if 'mean' in list(kwargs):
		ismean = kwargs['mean']
	else:
		ismean = 0
	if 'mask' in list(kwargs):
		mask = kwargs['mask']
	else:
		mask = np.ones((npts,),dtype=int)
	if 'parallel' in list(kwargs):
		isparallel = kwargs['paralllel']
	else:
		isparallel = 0
	if 'verbose' in list(kwargs):
		isverbose = kwargs['verbose']
	else:
		isverbose = 0

	# where is masking part?
	flags = (np.sum(mask[elements],1) == 3)

	# step1: get areas.
	Areas = np.array(GetAreas(elements[flags,:]+1, x, y));
	if ismean:
		Areas = Areas/np.sum(Areas)
	#weights = Areas/np.sum(Areas);
	#print(weights[weights < 0].shape)

	if isparallel:
		ncpus = multiprocessing.cpu_count()
		print('   initialize function.')
		temp_func = functools.partial(AreaSum,data=data,elements=elements,flags=flags, Areas=Areas)
		print('   initialize Pool.')
		print('   -- use tqdm')
		#with multiprocessing.Pool(8) as pool:
		#   output = tqdm.tqdm(pool.apply_async(temp_func,range(s[1])),
		#                           total=s[1])
		import p_tqdm
		output = p_tqdm.p_map(temp_func,range(s[1]))
	else:
		output = np.zeros((s[1],));
		for ii in tqdm.tqdm(range(s[1]),disable=~isverbose):
			output[ii] = np.nansum(Areas*np.mean(data[elements[flags,:],ii],axis=1))

	return output

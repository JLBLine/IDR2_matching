from numpy import *
from sys import path
path.append('/home/jline/Documents/cataloguing/PUMA/scripts')
import make_table_lib as mkl
import optparse
import copy
from astropy.table import Table, Column, MaskedColumn
from itertools import combinations

num_freqs = [20, 1, 1, 1, 1, 3]
pref_cats = ['nvss','sumss','comp_v11','atg20','vlssr','mrc']
matched_cats=['comp_v11','vlssr','mrc','sumss','nvss','atg20']
cat_freqs = [[76.0,84.0,92.0,99.0,107.0,115.0,122.0,130.0,143.0,151.0,158.0,166.0,174.0,181.0,189.0,197.0,204.0,212.0,220.0,227], [74.0], [408.0], [843.0],[1400.0], [5000.0,8000.0,20000.0]]
jstat_thresh = 0.1
chi_thresh = 10.0

def arcdist(RA1,RA2,Dec1,Dec2):  
	'''Calculates the arcdist between 2 points in deg'''
	dr = pi/180.0
	in1 = (90.0 - Dec1)*dr
	in2 = (90.0 - Dec2)*dr
	RA_d = (RA1 - RA2)*dr
	cosalpha = cos(in1)*cos(in2) + sin(in1)*sin(in2)*cos(RA_d)
	##Sometimes get floating point errors if two sources at exactly
	##the same position, so account for this:
	if cosalpha>1.0: cosalpha = 1.0
	elif cosalpha<-1.0: cosalpha = -1.0
	alpha = arccos(cosalpha)
	return alpha/dr


def make_the_table(t,sources,sources_stats):
	pos_round = 4
	flux_round = 4
	updated_ras = []
	updated_decs = []
	updated_rerrs = []
	updated_derrs = []

	for source in sources:
		##this is for choosing the most 'reliable' position as defined by
		##the input preference cats
		for p_cat in pref_cats:
			if p_cat in source.cats:
				pos_ind = source.cats.index(p_cat)
				break
		else:
			print 'no catalogue match in positional preferences - OMG'
		updated_ras.append(source.ras[pos_ind])
		updated_decs.append(source.decs[pos_ind])
		updated_rerrs.append(source.rerrs[pos_ind])
		updated_derrs.append(source.derrs[pos_ind])

	##Make an 'MWA' name based on position (limited in number of characters due to the RTS)
	##should change this really
	def make_name(ra,dec,ext):
		ra_str = mkl.deg_to_hour(ra)[:6]
		dec_str = mkl.deg_to_degmins(dec)[:5]
		if ext=='':
			return 'GLM J'+ra_str+dec_str
		else:
			return 'GLM J'+ra_str+dec_str+ext

	original_ras = [source.ras[0] for source in sources]
	original_decs = [source.decs[0] for source in sources]
	original_rerrs = [source.rerrs[0] for source in sources]
	original_derrs = [source.derrs[0] for source in sources]
	type_matches = [source.accept_type for source in sources_stats]
	
	##Name the sources based on position. If the sources were the result of a split,
	##name by the position of the base source, and add the letter that was defined by the
	##split letter in the stat report object
	names = []
	for i in xrange(len(type_matches)):
		if 'split' in type_matches[i]:
			ext = type_matches[i][-1]
			name = make_name(original_ras[i],original_decs[i],ext)
		else:
			name = make_name(original_ras[i],original_decs[i],'')
		names.append(name)

	##Create many columns of data
	t_name = Column(name='Name',data=names,description='Name based on IDR2 position of source',dtype=str)
	prim_names = [source.names[0] for source in sources]
	t_base_name = Column(name='IDR2_ID',data=prim_names,description='Name of %s component' %matched_cats[0],dtype=str)
	
	t_orra = Column(name='IDR2_RAJ2000',data=array(original_ras),description='IDR2 RA of source (deg)',unit='deg',dtype=float)
	t_ordec = Column(name='IDR2_DECJ2000',data=array(original_decs),description='IDR2 DEC of source (deg)',unit='deg',dtype=float)
	t_orrerr = Column(name='e_IDR2_RAJ2000',data=array(original_rerrs),description='IDR2 error on RA of source (deg)',unit='deg',dtype=float)
	t_orderr = Column(name='e_IDR2_DECJ2000',data=array(original_derrs),description='IDR2 error on DEC of source (deg)',unit='deg',dtype=float)
	##Add the columns
	t.add_columns([t_name,t_base_name,t_orra,t_ordec,t_orrerr,t_orderr])
	
	shift_ras = []
	shift_kgs_ras = []

	for i in xrange(len(updated_ras)):
		if updated_ras[i]> 180.0:
			shift_ras.append(updated_ras[i] -360.0)
		else:
			shift_ras.append(updated_ras[i])
		if original_ras[i]> 180.0:
			shift_kgs_ras.append(original_ras[i] -360.0)
		else:
			shift_kgs_ras.append(original_ras[i])

	cons_ra = sqrt(array(updated_rerrs)**2 + (array(shift_ras)-array(shift_kgs_ras))**2 + array(original_rerrs)**2)
	cons_dec = sqrt(array(updated_derrs)**2 + (array(updated_decs)-array(original_decs))**2 + array(original_derrs)**2)
	
	##Come up with some masks - for non-matches, mask out all PUMA related stuff, and the second type of error
	##For only o2 catalogues, mask out the SI errors
	num_cats = [len(set([cat for cat in source.cats if cat!='-100000.0'])) for source in sources]
	single_cat_mask = [num<=1 for num in num_cats]
	SI_err_mask = [num<=2 for num in num_cats]
	
	#t_e2err = MaskedColumn(name='e2_KGS_RAJ2000',data=array(cons_ra).round(pos_round),description='Conservative RA error estimate (deg)',unit='deg',mask=single_cat_mask,fill_value=None,dtype=float)
	#t_e2derr = MaskedColumn(name='e2_KGS_DECJ2000',data=array(cons_dec).round(pos_round),description='Conservative DEC error estimate (deg)',unit='deg',mask=single_cat_mask,fill_value=None,dtype=float)
	#t.add_columns([t_e2err,t_e2derr])
	
	##For just the KGS freq
	for cat in range(0,1):
		##See how many frequencies that source has
		num_freq = num_freqs[cat]
		##For every frequency, make a column of fluxes and flux errors, masking every value with -100000.0
		for freq in xrange(num_freq):
			fluxs = array([src.fluxs[cat][freq] for src in sources])
			ferrs = array([src.ferrs[cat][freq] for src in sources])
			t_flux = MaskedColumn(name='S_%d' %int(cat_freqs[cat][freq]),data=fluxs,description='Flux at %.1fMHz (Jy)' %float(cat_freqs[cat][freq]),mask=fluxs==-100000.0, fill_value=None,unit='Jy',dtype=float)
			t_ferr = MaskedColumn(name='e_S_%d' %int(cat_freqs[cat][freq]),data=ferrs,description='Flux error at %.1fMHz (Jy)' %float(cat_freqs[cat][freq]),mask=ferrs==-100000.0, fill_value=None,unit='Jy',dtype=float)
			t.add_columns([t_flux,t_ferr])
			
	#t_rel = Column(name='Reliability',data=ones(len(sources)),description='KATALOGSS Reliability classification',dtype=float)
	#t_mean = Column(name='Mean_Beam',data=ones(len(sources)),description='Average beam value of all detections at the source position (normalised)',dtype=float)
	#t_det = Column(name='N_det',data=ones(len(sources)),description='Number of detections made out of 71 observations',dtype=float)
	#t_exp = Column(name='N_exp',data=ones(len(sources)),description='Approximate number of detections expected out of 71 observations',dtype=float)
	#t_loc = Column(name='Local_Density',data=ones(len(sources)),description='Number of source candidates within 100 arcminutes',dtype=float)
	
	###Work out if v5 or other name to find what clustering was used
	#clust_flags = []
	
	#for name in prim_names:
		#if 'v5' in  name:
			#clust_flags.append('q')
		#else:
			#clust_flags.append('h')
			
	#t_clust = Column(name='R_cluster',data=array(clust_flags),description='The clustering radius use in the KGS process; q for one-quarter beam radius, h for one-half beam radius',dtype=str)
	
	#t.add_columns([t_rel,t_mean,t_det,t_exp,t_loc,t_clust])
	
	type_matches = [source.accept_type for source in sources_stats]
	#for i in xrange(len(type_matches)):
		##print type_matches[i]
		#if 'pos' in type_matches[i]:
			#type_matches[i] = 'isolated'
		#elif 'spec' in type_matches[i]:
			#type_matches[i] = 'dominant'
		#elif 'comb' in type_matches[i]:
			#type_matches[i] = 'multiple'
		#else:
			#pass
	t_stage = MaskedColumn(name='Match_Type',data=array(type_matches),description='The PUMA stage at which a decision was made',mask=single_cat_mask,fill_value=None,dtype=str)
	inspecs = [source.inspected for source in sources]
	inspecs = Column(name='Inspected',data=array(inspecs),description='Inspected flag. 0 if not visually inspected, 1 if catalogue and image information inspected, 2 if edited by user',dtype=int)
	
	t_numcats = MaskedColumn(name='Number_cats',data=array(num_cats),description='Number of catalogues in match',dtype=int)
	t_upra = MaskedColumn(name='Match_RAJ2000',data=array(updated_ras,dtype=float).round(pos_round),description='Updated RA of source (deg)',unit='deg',mask=single_cat_mask,fill_value=None,dtype=float)
	t_updec = MaskedColumn(name='Match_DECJ2000',data=array(updated_decs,dtype=float).round(pos_round),description='Updated DEC of source (deg)',unit='deg',mask=single_cat_mask,fill_value=None,dtype=float)
	t_uprerr = MaskedColumn(name='e_Match_RAJ2000',data=array(updated_rerrs,dtype=float).round(pos_round),description='Error on updated RA of source (deg)',unit='deg',mask=single_cat_mask,fill_value=None,dtype=float)
	t_upderr = MaskedColumn(name='e_Match_DECJ2000',data=array(updated_derrs,dtype=float).round(pos_round),description='Error on updated DEC of source (deg)',unit='deg',mask=single_cat_mask,fill_value=None,dtype=float)
	
	t.add_columns([t_stage,inspecs,t_numcats,t_upra,t_updec,t_uprerr,t_upderr])
	
	##create a mask for all of the errors reported for matched with only two cats - can't have an error on
	##a fit to just two data points, so statsmodel spits out nonsense
	SIs = array([source.SI for source in sources],dtype=float)
	SI_errs = array([source.SI_err for source in sources],dtype=float)
	t_SIs = MaskedColumn(name='SI',data=SIs.round(3),description='Spectral Index of Fit',mask=single_cat_mask,fill_value=None,dtype=float)
	t_SIerrs = MaskedColumn(name='e_SI',data=SI_errs.round(3),description='Std error on Spectral Index of Fit',mask=SI_err_mask,fill_value=None,dtype=float)
	t.add_columns([t_SIs,t_SIerrs])

	for cat in range(1,len(num_freqs)):
	##See how many frequencies that source has
		num_freq = num_freqs[cat]
		##For every frequency, make a column of fluxes and flux errors, masking every value with -100000.0
		for freq in xrange(num_freq):
			#for src in sources: print src.fluxs,src.fluxs[cat][freq]
			fluxs = array([src.fluxs[cat][freq] for src in sources],dtype=float)
			ferrs = array([src.ferrs[cat][freq] for src in sources],dtype=float)
			
			t_flux = MaskedColumn(name='S_%d' %int(cat_freqs[cat][freq]),data=fluxs.round(flux_round),description='Flux at %.1fMHz (Jy)' %float(cat_freqs[cat][freq]),mask=fluxs==-100000.0, fill_value=None,unit='Jy',dtype=float)
			t_ferr = MaskedColumn(name='e_S_%d' %int(cat_freqs[cat][freq]),data=ferrs.round(flux_round),description='Flux error at %.1fMHz (Jy)' %float(cat_freqs[cat][freq]),mask=ferrs==-100000.0, fill_value=None,unit='Jy',dtype=float)
			t.add_columns([t_flux,t_ferr])
	
	
	vlssr = [source.vlssr for source in sources]
	t_vlssr = MaskedColumn(name='VLSSr',data=array(vlssr),description='Names of all VLSSr matches, separated by commas',mask=vlssr=='', fill_value=None,dtype=str)
	mrc = [source.mrc for source in sources]
	t_mrc = MaskedColumn(name='MRC',data=array(mrc),description='Names of all MRC matches, separated by commas',mask=vlssr=='', fill_value=None,dtype=str)
	sumss = [source.sumss for source in sources]
	t_sumss = MaskedColumn(name='SUMSS',data=array(sumss),description='Names of all SUMSS matches, separated by commas',mask=vlssr=='', fill_value=None,dtype=str)
	nvss = [source.nvss for source in sources]
	t_nvss = MaskedColumn(name='NVSS',data=array(nvss),description='Names of all NVSS matches, separated by commas',mask=vlssr=='', fill_value=None,dtype=str)
	atg20 = [source.atg20 for source in sources]
	t_atg20 = MaskedColumn(name='ATG20',data=array(atg20),description='Names of all ATG20 matches, separated by commas',mask=vlssr=='', fill_value=None,dtype=str)
	t.add_columns([t_vlssr,t_mrc,t_sumss,t_nvss,t_atg20])
	
	accepts = [source.accept for source in sources]
	t_accepts = MaskedColumn(name='Outcomes',data=array(accepts),description='Accepted or rejected',mask=single_cat_mask,fill_value=None,dtype=str)
	chireds = [source.chi_resid for source in sources]
	t_chireds = MaskedColumn(name='Chi_sq_red',data=array(chireds),description='The chi_reduced value obtained from a wls fit',mask=SI_err_mask,fill_value=None,dtype=float)
	t.add_columns([t_accepts,t_chireds])

	return t

##Used to store source and group information
class source_group:
	def __init__(self):
		self.cats = []
		self.names = []
		self.ras = []
		self.rerrs = []
		self.decs = []
		self.derrs = []
		self.freqs = []
		self.fluxs = []
		self.ferrs = []
		self.majors = []
		self.minors = []
		self.PAs = []
		self.flags = []
		self.IDs = []
		self.SI = None
		self.intercept = None
		self.prob = None
		self.num_matches = None
		self.type_match = None
		self.SI_err = None
		self.intercept_err = None
		self.low_resids = None
		self.chi_resid = None
		self.jstat_resid = None
		self.vlssr = ''
		self.mrc = ''
		self.sumss = ''
		self.nvss = ''
		self.atg20 = ''
		self.inspected = -1
		self.accept = ''
		
def combine_flux(src_all,src_g,accepted_inds,num_matches):
	'''Takes a src_group() class that contains all group info. Indentifies which catalogue is repeated, combines
	the fluxes and fits a new line. Returns the reduced frequency list, the combined flux and flux error 
	arrays, as well as the line_fit object and residuals'''
	
	##Find repeated catalogues
	repeated_cats = set([src_all.cats[ind] for ind in accepted_inds if src_all.cats.count(src_all.cats[ind])>1])
	##This won't neccesarily be in the that the cats appear in src_all.cats so reorder
	repeat_indexs = [src_all.cats.index(cat) for cat in repeated_cats]
	repeated_cats = [cat for ind,cat in sorted(zip(repeat_indexs,repeated_cats),key=lambda pair: pair[0])]
	
	##These are used to test the combined spectrum
	temp_freqs = [src_all.freqs[i] for i in xrange(len(src_all.freqs)) if src_all.cats[i] not in repeated_cats]
	temp_fluxs = [src_all.fluxs[i] for i in xrange(len(src_all.fluxs)) if src_all.cats[i] not in repeated_cats]
	temp_ferrs = [src_all.ferrs[i] for i in xrange(len(src_all.ferrs)) if src_all.cats[i] not in repeated_cats]
	
	##Will need these for fitting/passing on to plotting function
	comb_freqs = []
	comb_fluxs = []
	comb_ferrs = []
	ra_ws = []
	dec_ws = []
	rerr_ws = []
	derr_ws = []
	
	##Need these in case of doing the split test
	resolved_diff_inds = []
	unrepeat_dists = []
	num_of_repeats = []
	
	for repeat_cat in repeated_cats:
		num_of_repeat = [i for i in xrange(len(src_all.names)) if (src_all.cats[i]==repeat_cat) and (i in accepted_inds)]
		num_of_repeats.append(num_of_repeat)
		
	##Need to put the source names that weren't repeated in their correct places
	for cat,name in zip(src_g.cats[1:],src_g.names[1:]):
		if cat not in repeated_cats and cat != '100000.0':
			if cat == 'vlssr':
				src_g.vlssr = name
			elif cat == 'mrc':
				src_g.mrc = name
			elif cat == 'sumss':
				src_g.sumss = name
			elif cat == 'nvss':
				src_g.nvss = name
			elif cat == 'atg20':
				src_g.atg20 = name
			else:
				pass
	
	##For each repeated catalogue:
	for repeat_cat in repeated_cats:
		##Find the frequency/ies of repeated cat
		comb_freq = src_all.freqs[src_all.cats.index(repeat_cat)]
		##Find the flux/es of the repeat_cat sources that were accepted by retained_sources()
		flux_to_comb = [src_all.fluxs[i] for i in xrange(len(src_all.fluxs)) if (src_all.cats[i]==repeat_cat) and (i in accepted_inds)]
		##Find the flux error/s of the repeat_cat sources that were accepted by retained_sources()
		ferr_to_comb = [src_all.ferrs[i] for i in xrange(len(src_all.ferrs)) if (src_all.cats[i]==repeat_cat) and (i in accepted_inds)]
		
		##ALso find all of the positional information to combine
		ras_to_comb = [src_all.ras[i] for i in xrange(len(src_all.ras)) if (src_all.cats[i]==repeat_cat) and (i in accepted_inds)]
		rerrs_to_comb = [src_all.rerrs[i] for i in xrange(len(src_all.rerrs)) if (src_all.cats[i]==repeat_cat) and (i in accepted_inds)]
		decs_to_comb = [src_all.decs[i] for i in xrange(len(src_all.decs)) if (src_all.cats[i]==repeat_cat) and (i in accepted_inds)]
		derrs_to_comb = [src_all.derrs[i] for i in xrange(len(src_all.derrs)) if (src_all.cats[i]==repeat_cat) and (i in accepted_inds)]
		
		##This is to write down all the names of the sources combined
		names_to_comb = [src_all.names[i] for i in xrange(len(src_all.names)) if (src_all.cats[i]==repeat_cat) and (i in accepted_inds)]
		name_string = ''
		for name in names_to_comb: name_string += ','+name
		
		if repeat_cat == 'vlssr':
			src_g.vlssr = name_string[1:]
		elif repeat_cat == 'mrc':
			src_g.mrc = name_string[1:]
		elif repeat_cat == 'sumss':
			src_g.sumss = name_string[1:]
		elif repeat_cat == 'nvss':
			src_g.nvss = name_string[1:]
		elif repeat_cat == 'atg20':
			src_g.atg20 = name_string[1:]
		else:
			pass
		
		##TEST TO SEE IF THE REPEATED SOURCES ARE RESOLVED (BY A GIVEN RESOLUTION THRESHOLD)
		##Need to do the split test here even if not propagating to the final catalogue
		##-------------------------------------------------------------------------------------------------------
		
		##Give it a value even if not splitting
		#if split==0:
			#dist_test = 0.020833333 ##1.25 arcmin
		#else:
			#dist_test = dec_to_deg(split)
		#big_inds = []
		
		#n = len(ras_to_comb)
		
		#for i in range(0,n+1):
			#for j in range(i+1,n):
				#dist = arcdist(ras_to_comb[i],ras_to_comb[j],decs_to_comb[i],decs_to_comb[j])
				#if dist>dist_test:
					#big_inds.append([i,j])
		#resolved_diff_inds.append(big_inds)
		##-------------------------------------------------------------------------------------------------------
		
		comb_flux = sum(flux_to_comb) #Sum the fluxes
		comb_ferr = zeros(1)
		for ferr in ferr_to_comb: comb_ferr+= ferr**2 #Add the errors in quadrature
		comb_ferr = comb_ferr**0.5

		##Need these later to fit the combined flux and to populate a new src_g if
		##the fit passes
		
		temp_freqs.append(comb_freq)
		temp_fluxs.append([comb_flux])
		temp_ferrs.append([comb_ferr[0]])
		comb_freqs.append(comb_freq)
		comb_fluxs.append(comb_flux)
		comb_ferrs.append(comb_ferr)
		
		##Weight by the first flux in the flux list (in case catalogue has multiple frequencies)
		flux_s = [flux[0] for flux in flux_to_comb]

		##A a little code in case sources are very close to RA=0, and some are reporting 359.9,
		##and others 0.001 - if you don't account for this, the weigthed position goes mental
		wrap = 'no'
		for combo in combinations(ras_to_comb,2):
			diff = combo[0]-combo[1]
			if abs(diff)>180.0:
				wrap='yes'
		if wrap=='yes':
			for i in xrange(len(ras_to_comb)):
				if ras_to_comb[i]<180.0:
					ras_to_comb[i]+=360.0
		##Weight the sources by their flux
		weights = [flux/sum(flux_s) for flux in flux_s]

		##Do the weighting, if weighted RA is above 360 deg, rescale, 
		##add the errors as shown in the write up
		ra_w = dot(ras_to_comb,weights)
		if ra_w>360.0: ra_w-=360.0
		dec_w = dot(decs_to_comb,weights)
		rerr_w = (dot(rerrs_to_comb,weights)**2)**0.5
		derr_w = (dot(derrs_to_comb,weights)**2)**0.5
		
		ra_ws.append(ra_w)
		dec_ws.append(dec_w)
		rerr_ws.append(rerr_w)
		derr_ws.append(derr_w)

	log_temp_freqs = []
	log_temp_fluxs = []
	log_temp_ferrs = []
	
	##Get the sources out of the array in list format (which is used later when making the sources
	##to add to the final table)
	for i in xrange(len(temp_freqs)):
		for j in xrange(len(temp_freqs[i])):
			if temp_fluxs[i][j] == -100000.0 or isnan(temp_fluxs[i][j])==True:
				pass
			else:
				log_temp_freqs.append(log(temp_freqs[i][j]))
	for i in xrange(len(temp_freqs)):
		for j in xrange(len(temp_freqs[i])):
			if temp_fluxs[i][j] == -100000.0 or isnan(temp_fluxs[i][j])==True:
				pass
			else:
				log_temp_fluxs.append(log(temp_fluxs[i][j]))
	for i in xrange(len(temp_freqs)):
		for j in xrange(len(temp_freqs[i])):
			if temp_fluxs[i][j] == -100000.0 or isnan(temp_fluxs[i][j])==True:
				pass
			else:
				log_temp_ferrs.append(temp_ferrs[i][j]/temp_fluxs[i][j])

	##Fit and find residuals to the combined spectrum
	comb_fit,comb_jstat,comb_bse,comb_chi_red = mkl.fit_line(array(log_temp_freqs),array(log_temp_fluxs),array(log_temp_ferrs))
	
	##Find out where in srg_g the repeated cats appear
	repeat_cat_inds = [src_g.cats.index(cat) for cat in repeated_cats]
	
	split_flag=''
	##Make labels for when we're plotting a put in combined_names
	combined_names = []
	if comb_jstat<=jstat_thresh or comb_chi_red<=chi_thresh:
		##Create the combined source no matter what for plotting purposes
		##Loop over all the combined sources and repopulate the entries of a src_g
		##at the point where the repeated catalogues appear
		for i in xrange(len(comb_fluxs)):
			srcg_ind = repeat_cat_inds[i]
			src_g.ras[srcg_ind] = ra_ws[i]
			src_g.rerrs[srcg_ind] = rerr_ws[i]
			src_g.decs[srcg_ind] = dec_ws[i]
			src_g.derrs[srcg_ind] = derr_ws[i]
			src_g.PAs[srcg_ind] = -100000.0
			src_g.majors[srcg_ind] = -100000.0
			src_g.minors[srcg_ind] = -100000.0
			src_g.names[srcg_ind] = "Combined-%s" %src_g.cats[srcg_ind]
			combined_names.append("Combined-%s" %src_g.cats[srcg_ind])
			src_g.fluxs[srcg_ind] = [comb_fluxs[i]]
			src_g.ferrs[srcg_ind] = [comb_ferrs[i]]
		#srg_g.freqs = temp_freqs
		src_g.SI = float(comb_fit.params[0])
		src_g.intercept = comb_fit.params[1]
		src_g.SI_err = comb_bse[0]
		src_g.intercept_err = comb_bse[1]
		src_g.chi_resid = comb_chi_red
		src_g.epsilon_red = comb_jstat
	
		##If good fit, report that in the final stats object
		if comb_chi_red<=2:
			src_g.low_resids = 0
		else:
			src_g.low_resids = 1
			
		dom_crit = 'Accepted -\ncombined'
		return dom_crit, [src_g], comb_jstat, comb_chi_red

	##If it fails, still return all the info to the plot so we can see what's going on
	else:
		dom_crit = 'To eyeball'
		return dom_crit, 'nyope', comb_jstat, comb_chi_red
		
def get_srcg(info):
	'''Takes a string which contains the information for all sources in particular combination
	Uses num_freqs to work out where each piece of information is, and then return the relevant
	information in a source_group() class. Any catalogues with no matches are entered as -10000.0'''
	src_g = source_group()
	##Work out where in the string each different catalogue will start, using num_freqs
	##to work out how many entries each catalogue will have
	indexes = [(14+((i-1)*3)) for i in num_freqs]
	starts = [0]
	for i in xrange(len(indexes)-1): starts.append(sum(indexes[:i+1]))
	#all_freqs = []
	#all_fluxs = []
	#all_ferrs = []
	for j in xrange(len(starts)): 
		num_freq = num_freqs[j]
		ind = starts[j]
		#cat = info[ind]
		freqss = []
		fluxss = []
		ferrss = []
		for k in xrange(num_freq):
			##Test to see if flux is a nan or -100000.0; make sure all flux/freq info is -100000.0 if so
			if isnan(float(info[7+ind+(3*k)])) == False or float(info[7+(3*i)]) != -100000.0:
				freqss.append(float(info[6+ind+(3*k)]))
				fluxss.append(float(info[7+ind+(3*k)]))
				ferrss.append(float(info[8+ind+(3*k)]))
			else:
				freqss.append(-100000.0)
				fluxss.append(-100000.0)
				ferrss.append(-100000.0)
			#if float(info[6+ind+(3*k)])!=-100000.0: all_freqs.append(float(info[6+ind+(3*k)]))
			#if float(info[7+ind+(3*k)])!=-100000.0: all_fluxs.append(float(info[7+ind+(3*k)]))
			#if float(info[8+ind+(3*k)])!=-100000.0: all_ferrs.append(float(info[8+ind+(3*k)]))
		src_g.freqs.append(freqss)
		src_g.fluxs.append(fluxss)
		src_g.ferrs.append(ferrss)
		src_g.cats.append(info[ind])
		src_g.names.append(info[ind+1])
		src_g.ras.append(float(info[ind+2]))
		src_g.rerrs.append(float(info[ind+3]))
		src_g.decs.append(float(info[ind+4]))
		src_g.derrs.append(float(info[ind+5]))
		src_g.majors.append(float(info[ind+9+((num_freq-1)*3)]))
		src_g.minors.append(float(info[ind+10+((num_freq-1)*3)]))
		src_g.PAs.append(float(info[ind+11+((num_freq-1)*3)]))
		src_g.flags.append(info[ind+12+((num_freq-1)*3)])
		src_g.IDs.append(info[ind+13+((num_freq-1)*3)])
		src_g.prob = float(info[-1])
	return src_g

def get_allinfo(all_info):
	'''Takes a list of strings. Each string is a line containing the information for a single source
	in a matched group from an output file of calculate_bayes.py. Gets all of the information
	from each entry and returns them in a source_group() class'''
	src_all = source_group()
	for entry in all_info:
		info = entry.split()
		src_all.cats.append(info[0])
		src_all.names.append(info[1])
		src_all.ras.append(float(info[2]))
		src_all.rerrs.append(float(info[3]))
		src_all.decs.append(float(info[4]))
		src_all.derrs.append(float(info[5]))
		src_all.majors.append(float(info[-5]))
		src_all.minors.append(float(info[-4]))
		src_all.PAs.append(float(info[-3]))
		src_all.flags.append(info[-2])
		src_all.IDs.append(info[-1])
		
		##If the source only has one frequency. Append as an array
		##so that all entries to src_all.freqs etc are of the same
		##type. This deals with cats with multiple freqs, otherwise
		##the position, name etc will have to be repeated for each
		##frequency
		if len(info)==14:
			src_all.freqs.append(array([float(info[6])]))          
			src_all.fluxs.append(array([float(info[7])]))
			src_all.ferrs.append(array([float(info[8])]))
			#src_all.freqs.append(float(info[6]))  ##Left here in case ever want to append just the freq, not
			#src_all.fluxs.append(float(info[7]))  ##an array
			#src_all.ferrs.append(float(info[8]))
			
		##If not, work out how many freqs there are and append to lists
		else:
			extra = (len(info)-14) / 3
			freqs = []
			fluxs = []
			ferrs = []
			for i in xrange(extra+1):
				##Test to see if flux is a nan or -100000.0; make sure all flux/freq info is -100000.0 if so
				if isnan(float(info[7+(3*i)])) == False or float(info[7+(3*i)]) != -100000.0:
					freqs.append(float(info[6+(3*i)]))
					fluxs.append(float(info[7+(3*i)]))
					ferrs.append(float(info[8+(3*i)]))
				else:
					freqs.append(-100000.0)
					fluxs.append(-100000.0)
					ferrs.append(-100000.0)
			src_all.freqs.append(array(freqs))
			src_all.fluxs.append(array(fluxs))
			src_all.ferrs.append(array(ferrs))
				#src_all.freqs.append(float(info[6+(3*i)]))
				#src_all.fluxs.append(float(info[7+(3*i)]))
				#src_all.ferrs.append(float(info[8+(3*i)]))
	return src_all
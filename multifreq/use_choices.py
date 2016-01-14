#!/usr/bin/python
import numpy as np
from itertools import combinations
import optparse
from sys import path
#path.append('/home/jline/Documents/cataloguing/PUMA/scripts')
import make_table_lib as mkl
import subprocess
from astropy.table import Table, Column, MaskedColumn
from make_table_function import *

dr= np.pi/180.0

parser = optparse.OptionParser()

parser.add_option('-f', '--cat_freqs', 
	help='Enter the frequencies of each catalogue')

parser.add_option('-p', '--pref_cats',
	help='Enter names of matched cataloges, in order of preferable position')

parser.add_option('-m', '--matched_cats', 
	help='Enter the names of each catalogue, base catalogue first')

parser.add_option('-i', '--input_bayes', 
	help='Enter name of eyeball bayes file')

parser.add_option('-t', '--tag', 
	help='Enter tag for output files')

parser.add_option('-d', '--decisions', 
	help='Enter names and decisions of the user')

parser.add_option('-a', '--accept', action='store_true', default=False,
	help='Enter name for output catalogue')

options, args = parser.parse_args()

cat_fs = options.cat_freqs.split(',')
cat_freqs= []
for fs in cat_fs:
	split = fs.split('~')
	split = map(float,split)
	cat_freqs.append(split)

matched_cats = options.matched_cats.split(',')

pref_cats = []
for pref in options.pref_cats.split(','): pref_cats.append(pref)

num_freqs = []
for freq in cat_fs: num_freqs.append(len(freq.split('~')))

closeness = 1.17/60.0
high_prob = 0.95
low_prob = 0.8
chi_thresh = 10.0
jstat_thresh = 0.1
num_freqs = num_freqs
split = 0

mkl.closeness = 1.17/60.0
mkl.high_prob = 0.95
mkl.low_prob = 0.8
mkl.chi_thresh = 10.0
mkl.jstat_thresh = 0.1
mkl.num_freqs = num_freqs
mkl.split = 0

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def deg_to_degmins(x,style):	    #converts angle degrees form in to dd:mm:ss.ss
	x=float(x)
	deg=abs(x)
	degr=deg-int(deg)
	mins=(degr)*60.00
	secs=(mins-int(mins))*60.0
	if mins!=0:
		if -1e-5<=(secs-60)<1e-5:
			mins=mins+1
			secs=0.0
	if style == 'info':
		if x>0:
			return '+%02d %02d %04.1f' %(int(deg),int(mins),secs)
		if x<0:
			return '-%02d %02d %04.1f' %(int(deg),int(mins),secs)
	elif style == 'name':
		if x>0:
			return '+%02d%02d%04.1f' %(int(deg),int(mins),secs)
		if x<0:
			return '-%02d%02d%04.1f' %(int(deg),int(mins),secs)

def deg_to_hour(x,style):    #converts angle in degrees in to hh:mm:ss.ss, must input as a string
	x=float(x)
	deg=abs(x)
	hr=deg/15.0
	mins=(hr-int(hr))*60.0
	secs=(mins-int(mins))*60.0
	if mins!=0:
		if -1e-5<=(secs-60)<1e-5:
			mins=mins+1
			secs=0.0
	if style == 'info':
		if x>0:
			return '%02d %02d %04.1f' %(int(hr),int(mins),secs)
		if x<0:
			return '-%02d %02d %04.1f' %(int(hr),int(mins),secs)
	elif style == 'name':
		if x>0:
			return '%02d%02d%04.1f' %(int(hr),int(mins),secs)
		if x<0:
			return '-%02d%02d%04.1f' %(int(hr),int(mins),secs)
	
	
def xtick_RA(x):    #converts angle in degrees in to hh:mm:ss.ss, must input as a string
	x=float(x)
	deg=abs(x)
	hr=deg/15.0
	mins=(hr-int(hr))*60.0
	secs=(mins-int(mins))*60.0
	if mins!=0:
		if -1e-5<=(secs-60)<1e-5:
			mins=mins+1
			secs=0.0
	if x>0:
		return '%02d:%02d:%02d' %(int(hr),int(mins),secs)
	if x<0:
		return '-%02d:%02d:%02d' %(int(hr),int(mins),secs)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def combine_all_flux(src_all,src_g):
	'''Takes a src_group() class that contains all group info. Indentifies which catalogue is repeated, combines
	the fluxes and fits a new line. Returns the reduced frequency list, the combined flux and flux error 
	arrays, as well as the line_fit object and residuals'''
	
	accepted_inds = xrange(len(src_all.names))
	
	##Find repeated catalogues
	repeated_cats = set([src_all.cats[ind] for ind in accepted_inds if src_all.cats.count(src_all.cats[ind])>1])
	##This won't neccesarily be in the that the cats appear in src_all.cats so reorder
	repeat_indexs = [src_all.cats.index(cat) for cat in repeated_cats]
	repeated_cats = [cat for ind,cat in sorted(zip(repeat_indexs,repeated_cats),key=lambda pair: pair[0])]
	
	##These are used to test the combined spectrum
	temp_freqs = [src_all.freqs[i] for i in xrange(len(src_all.freqs)) if src_all.cats[i] not in repeated_cats]
	temp_fluxs = [src_all.fluxs[i] for i in xrange(len(src_all.fluxs)) if src_all.cats[i] not in repeated_cats]
	temp_ferrs = [src_all.ferrs[i] for i in xrange(len(src_all.ferrs)) if src_all.cats[i] not in repeated_cats]
	
	base_freqs = [src_all.freqs[i] for i in xrange(len(src_all.freqs)) if src_all.cats[i] not in repeated_cats]
	base_fluxs = [src_all.fluxs[i] for i in xrange(len(src_all.fluxs)) if src_all.cats[i] not in repeated_cats]
	base_ferrs = [src_all.ferrs[i] for i in xrange(len(src_all.ferrs)) if src_all.cats[i] not in repeated_cats]
	base_cats = [src_all.cats[i] for i in xrange(len(src_all.ferrs)) if src_all.cats[i] not in repeated_cats]
	base_ras = [src_all.ras[i] for i in xrange(len(src_all.ferrs)) if src_all.cats[i] not in repeated_cats]
	base_rerrs = [src_all.rerrs[i] for i in xrange(len(src_all.ferrs)) if src_all.cats[i] not in repeated_cats]
	base_decs = [src_all.decs[i] for i in xrange(len(src_all.ferrs)) if src_all.cats[i] not in repeated_cats]
	base_derrs = [src_all.derrs[i] for i in xrange(len(src_all.ferrs)) if src_all.cats[i] not in repeated_cats]
	base_majors = [src_all.majors[i] for i in xrange(len(src_all.ferrs)) if src_all.cats[i] not in repeated_cats]
	base_minors = [src_all.minors[i] for i in xrange(len(src_all.ferrs)) if src_all.cats[i] not in repeated_cats]
	base_PAs = [src_all.PAs[i] for i in xrange(len(src_all.ferrs)) if src_all.cats[i] not in repeated_cats]
	base_names = [src_all.names[i] for i in xrange(len(src_all.ferrs)) if src_all.cats[i] not in repeated_cats]
	
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
		
		comb_flux = sum(flux_to_comb) #Sum the fluxes
		comb_ferr = np.zeros(1)
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
		ra_w = np.dot(ras_to_comb,weights)
		if ra_w>360.0: ra_w-=360.0
		dec_w = np.dot(decs_to_comb,weights)
		rerr_w = (np.dot(rerrs_to_comb,weights)**2)**0.5
		derr_w = (np.dot(derrs_to_comb,weights)**2)**0.5
		
		ra_ws.append(ra_w)
		dec_ws.append(dec_w)
		rerr_ws.append(rerr_w)
		derr_ws.append(derr_w)
	
	log_temp_freqs = []
	log_temp_fluxs = []
	log_temp_ferrs = []
	
	##Get the sources out of the array in list format (which is used later when making the sources
	##to add to the final table)
	for freqs in temp_freqs:
		for freq in freqs: log_temp_freqs.append(np.log(freq))
	for fluxs in temp_fluxs:
		for flux in fluxs: log_temp_fluxs.append(np.log(flux))
	for i in xrange(len(temp_ferrs)):
		for j in xrange(len(temp_ferrs[i])): log_temp_ferrs.append(temp_ferrs[i][j]/temp_fluxs[i][j])

	##Fit and find residuals to the combined spectrum
	comb_fit,comb_jstat,comb_bse,comb_chi_red = mkl.fit_line(np.array(log_temp_freqs),np.array(log_temp_fluxs),np.array(log_temp_ferrs))
	
	##Find out where in srg_g the repeated cats appear
	repeat_cat_inds = [src_g.cats.index(cat) for cat in repeated_cats]
	
	split_flag=''
	##Make labels for when we're plotting a put in combined_names
	combined_names = []
	#if comb_jstat<=jstat_thresh or comb_chi_red<=chi_thresh:
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
		src_g.ferrs[srcg_ind] = comb_ferrs[i]
		
	##Need to make sure the single sources in src_all are used rather than
	##whatever was in src_g
	for i in xrange(len(base_fluxs)):
		srcg_ind = src_g.cats.index(base_cats[i])
		src_g.ras[srcg_ind] = base_ras[i]
		src_g.rerrs[srcg_ind] = base_rerrs[i]
		src_g.decs[srcg_ind] = base_decs[i]
		src_g.derrs[srcg_ind] = base_derrs[i]
		src_g.PAs[srcg_ind] = base_PAs[i]
		src_g.majors[srcg_ind] = base_majors[i]
		src_g.minors[srcg_ind] = base_minors[i]
		src_g.names[srcg_ind] = base_names[i]
		src_g.fluxs[srcg_ind] = base_fluxs[i]
		src_g.ferrs[srcg_ind] = base_ferrs[i]
		
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
			
	return src_g

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
		dist_test = 0.020833333 ##1.25 arcmin
		big_inds = []
		
		n = len(ras_to_comb)
		
		for i in range(0,n+1):
			for j in range(i+1,n):
				dist = mkl.arcdist(ras_to_comb[i],ras_to_comb[j],decs_to_comb[i],decs_to_comb[j])
				if dist>dist_test:
					big_inds.append([i,j])
		resolved_diff_inds.append(big_inds)
		##-------------------------------------------------------------------------------------------------------
		
		comb_flux = sum(flux_to_comb) #Sum the fluxes
		comb_ferr = np.zeros(1)
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
		ra_w = np.dot(ras_to_comb,weights)
		if ra_w>360.0: ra_w-=360.0
		dec_w = np.dot(decs_to_comb,weights)
		rerr_w = (np.dot(rerrs_to_comb,weights)**2)**0.5
		derr_w = (np.dot(derrs_to_comb,weights)**2)**0.5
		
		ra_ws.append(ra_w)
		dec_ws.append(dec_w)
		rerr_ws.append(rerr_w)
		derr_ws.append(derr_w)
	
	
	##Flag distance between repeated cats as being larger than designated resolution
	big_flags = [0 for i in xrange(len(resolved_diff_inds))]
	for i in xrange(len(resolved_diff_inds)):
		if len(resolved_diff_inds[i])>0: big_flags[i] = 1 
	##Now each repeat_cat has a 1 flag if large separation, 0 if not
	
	set_freqs = []
	set_fluxs = []
	set_ferrs = []
	set_fits = []
	set_jstat = []
	set_bse = []
	set_red = []
	set_cats = []
	set_names = []
	big_sep = 'no'
	##If all repeated cats have large separation, and they have the same amount of repeated sources:
	if 0 not in big_flags and len(list(set([len(reap) for reap in num_of_repeats])))==1:
		##If more than one repeated cat (need more than one data point to get some spectral info)
		if len(repeated_cats)>1:
			##Find the 'sub' set matches, so sources that could be combined to make components
			sets = []
			#name_sets = []
			for src in num_of_repeats[0]:
				match = [src]
				#names = [src_all.names[src]]
				for other_srcs in num_of_repeats[1:]:
					for other_src in other_srcs:
						if mkl.arcdist(src_all.ras[src],src_all.ras[other_src],src_all.decs[src],src_all.decs[other_src])<dist_test:
							match.append(other_src)
							#names.append(src_all.names[other_src])
				sets.append(match)
				#name_sets.append(names)
			##If all the sets found have the same amount of components, and only have one source from
			##each repeated catalogue
			if len(list(set([len(sset) for sset in sets])))==1 and len(sets[0])==len(repeated_cats):
				big_sep = 'yes'
				for sset in sets:
					freqs = [src_all.freqs[src][0] for src in sset]
					fluxs = [src_all.fluxs[src][0] for src in sset]
					ferrs = [src_all.ferrs[src][0] for src in sset]
					names = [src_all.names[src] for src in sset]
					cats = [src_all.cats[src] for src in sset]
					set_freqs.append(freqs)
					set_fluxs.append(fluxs)
					set_ferrs.append(ferrs)
					set_names.append(names)
					set_cats.append(cats)
				
				flux_to_weight = [src_all.fluxs[i][0] for i in xrange(len(src_all.fluxs)) if (src_all.cats[i] not in repeated_cats)]
				freq_to_weight = [src_all.freqs[i][0] for i in xrange(len(src_all.freqs)) if (src_all.cats[i] not in repeated_cats)]
				ferr_to_weight = [src_all.ferrs[i][0] for i in xrange(len(src_all.ferrs)) if (src_all.cats[i] not in repeated_cats)]
				cats_to_weight = [src_all.cats[i] for i in xrange(len(src_all.cats)) if (src_all.cats[i] not in repeated_cats)]
				names_to_weight = [src_all.names[i] for i in xrange(len(src_all.names)) if (src_all.cats[i] not in repeated_cats)]
				
				##Find all the fluxs of the repeated cats, come up with weights for the single sources based on each
				##individual repeated catalogue, then take an average of these weights
				fluxs_for_weights = [[src_all.fluxs[sset[src]][0] for sset in sets] for src in xrange(len(sets[0]))]
				fluxs_weights = [[flux/sum(fluxs) for flux in fluxs] for fluxs in fluxs_for_weights]
				flux_weights = np.array([np.mean([[weights[weight]] for weights in fluxs_weights]) for weight in xrange(len(fluxs_weights[0]))])
				
				##For each set of freq,fluxs in the new set matched, append the weighted freq, flux and ferr of the
				##sources that have been split up
				for i in xrange(len(set_freqs)):
					weighted_fluxs = np.array(flux_to_weight)*flux_weights[i]
					weighted_errs = np.array(ferr_to_weight)*flux_weights[i]
					for j in xrange(len(weighted_fluxs)):
						set_freqs[i].append(freq_to_weight[j])
						set_fluxs[i].append(weighted_fluxs[j])
						set_ferrs[i].append(weighted_errs[j])
						set_cats[i].append(cats_to_weight[j])
						set_names[i].append(names_to_weight[j])

				##For every set of frequencies, ferrs, names and fluxes in the set, order the fluxes, names and ferrs by the frequencies
				set_fluxs = [[flux for flux,freq in sorted(zip(fluxs,freqs), key=lambda pair: pair[1]) ] for fluxs,freqs in zip(set_fluxs,set_freqs)]
				set_ferrs = [[ferr for ferr,freq in sorted(zip(ferrs,freqs), key=lambda pair: pair[1]) ] for ferrs,freqs in zip(set_ferrs,set_freqs)]
				set_cats = [[cat for cat,freq in sorted(zip(cats,freqs), key=lambda pair: pair[1]) ] for cats,freqs in zip(set_cats,set_freqs)]
				set_names = [[name for name,freq in sorted(zip(names,freqs), key=lambda pair: pair[1]) ] for names,freqs in zip(set_names,set_freqs)]
				set_freqs = [sorted(freq) for freq in set_freqs]

			for i in xrange(len(set_fluxs)):
				freqs = set_freqs[i]
				fluxs = set_fluxs[i]
				ferrs = set_ferrs[i]
				fit,jstat,bse,red = mkl.fit_line(np.log(freqs),np.log(fluxs),np.array(ferrs)/np.array(fluxs))
				set_fits.append(fit)
				set_jstat.append(jstat)
				set_bse.append(bse)
				set_red.append(red)
	log_temp_freqs = []
	log_temp_fluxs = []
	log_temp_ferrs = []
	
	##Get the sources out of the array in list format (which is used later when making the sources
	##to add to the final table)
	for freqs in temp_freqs:
		for freq in freqs: log_temp_freqs.append(np.log(freq))
	for fluxs in temp_fluxs:
		for flux in fluxs: log_temp_fluxs.append(np.log(flux))
	for i in xrange(len(temp_ferrs)):
		for j in xrange(len(temp_ferrs[i])): log_temp_ferrs.append(temp_ferrs[i][j]/temp_fluxs[i][j])

	##Fit and find residuals to the combined spectrum
	comb_fit,comb_jstat,comb_bse,comb_chi_red = mkl.fit_line(np.array(log_temp_freqs),np.array(log_temp_fluxs),np.array(log_temp_ferrs))
	
	##Find out where in srg_g the repeated cats appear
	repeat_cat_inds = [src_g.cats.index(cat) for cat in repeated_cats]
	
	split_flag=''
	##Make labels for when we're plotting a put in combined_names
	combined_names = []
	#if comb_jstat<=jstat_thresh or comb_chi_red<=chi_thresh:
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
		src_g.ferrs[srcg_ind] = comb_ferrs[i]
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
	
	return [src_g] 

reject_pile = open(options.tag+'_reject_pile.txt','w+')
reject_list = []
reject_list_file = open(options.tag+'_reject_list.txt','w+')
#v5_pile = open(options.tag+'_v5_pile.txt','w+')
#v5_list = []
#v5_list_file = open(options.tag+'_v5_list.txt','w+')
sources = []
sources_stats = []

##Input eyeball file
bayes_comp = open(options.input_bayes).read().split('END_GROUP')
del bayes_comp[-1]

num_found = 0
print len(bayes_comp)

for comp in bayes_comp:

	##Get the information into nice usable forms, and get rid of empty/pointless
	##entries
	chunks = comp.split('START_COMP')
	all_info = chunks[0].split('\n')
	
	##FOR SOME REASON CAN'T DO BOTH OF THESE LINES IN THE SAME FOR LOOP?!?!?!
	for entry in all_info:   
		if entry=='': del all_info[all_info.index(entry)]
	for entry in all_info:
		if 'START' in entry: del all_info[all_info.index(entry)]
		
	matches = chunks[1].split('\n')
	
	del matches[0],matches[-1]
	stats = matches[-1]
	del matches[-2:]
	
	match1 = matches[0].split()
	src_g = get_srcg(match1)
	
	##Get some info and find which catalogues are present 
	src_all = get_allinfo(all_info)
	
	answers = open(options.decisions,'r').read().split('\n')
	if answers[-1] == '': del answers[-1]
	edit_names = [info.split()[0] for info in answers]
	
	#print len(edit_names)
	
	#print edit_names
	
	if src_all.names[0] in edit_names:
		for info in answers:
			if src_all.names[0] == info.split()[0]:
				num_found += 1
				answer = info.split()[1]
				try: answer_2 = info.split()[2]
				except: answer_2 = 'meh'
				
		#print src_all.names,answer
				
		meh,num_matches,accept_matches,accepted_inds,accept_type,stage = stats.split()
		g_stats = mkl.group_stats()
		g_stats.num_matches = int(num_matches)
		g_stats.retained_matches = int(accept_matches)
		g_stats.accept_type = stage
		g_stats.outcome = accept_type
		##See how many unique catalogues there are
		num_cats = len(set([cat for cat in src_all.cats if cat!='-100000.0']))
		g_stats.num_cats = num_cats
		
		num_matches = len(matches)
		if answer == 'r':
			if comp[0] == 'S':
				pass
			else:
				comp = comp[1:]
			reject_list.append(src_all.names[0])
			reject_pile.write(comp)
			reject_pile.write('END_GROUP\n')
		##This is if we disagree with PUMA and want to combine every matched source
		elif answer == 'c':

			present_cats = [cat for cat in src_all.cats if cat!='-100000.0']
	
			flags = src_all.flags[0].split('+')
				
			meh,num_matches,accept_matches,accepted_inds,accept_type,stage = stats.split()
			g_stats = mkl.group_stats()
			g_stats.num_matches = int(num_matches)
			g_stats.retained_matches = int(accept_matches)
		
			src_comb = combine_all_flux(src_all,src_g)
				
			num_matches = len(matches)
			src_comb.inspected = 2
			src_comb.accept = accept_type
			src_comb.accept_type = stage
			sources.append(src_comb)
			#g_stats.accept_type = 'multiple'
		
			g_stats.outcome = accept_type
			if 'combine' in stage: 
				stage = 'multiple'
			elif 'position' in stage:
				stage = 'isolated'
			g_stats.accept_type = stage
	
			##See how many unique catalogues there are
			num_cats = len(set([cat for cat in src_all.cats if cat!='-100000.0']))
			g_stats.num_cats = num_cats

			sources_stats.append(g_stats)

			##Need this bit for plotting later on
			if comp[0] == 'S':
				pass
			else:
				comp = comp[1:]
			reject_list.append(src_all.names[0])
			reject_pile.write(comp)
			reject_pile.write('END_GROUP\n')
	
			
		elif answer == 'a':
			if answer_2 == 'c':
				comb_sources = combine_flux(src_all,src_g,map(int,accepted_inds.split(',')),num_matches)
				comb_sources[0].inspected = 1
				comb_sources[0].accept = accept_type
				sources.append(comb_sources[0])
				
				#g_stats.accept_type = 'check-combine'
				sources_stats.append(g_stats)
				#break
			else:
				accepted_match = matches[int(answer_2)-1].split()
				jstat_resids,params,bses,chi_resids = mkl.calculate_resids([accepted_match])
				source = get_srcg(accepted_match)
				source.SI = params[0][0]
				source.intercept = params[0][1]
				source.SI_err = bses[0][0]
				source.intercept_err = bses[0][1]
				source.chi_resid = chi_resids[0]
				source.epsilon_red = jstat_resids[0]
				source.inspected = 1
				if chi_resids[0]<=2:
					source.low_resids = 0
				else:
					source.low_resids = 1
					
				for cat,name in zip(source.cats[1:],source.names[1:]):
					if cat == 'vlssr':
						source.vlssr = name
					elif cat == 'mrc':
						source.mrc = name
					elif cat == 'sumss':
						source.sumss = name
					elif cat == 'nvss':
						source.nvss = name
					elif cat == 'atg20':
						source.atg20 = name
					else:
						pass
				source.accept = accept_type
				sources.append(source)
				#g_stats.accept_type = 'check-single'
				sources_stats.append(g_stats)
		elif answer == 'v5':
			if comp[0] == 'S':
				pass
			else:
				comp = comp[1:]
			v5_list.append(src_all.names[0])
			v5_pile.write(comp)
			v5_pile.write('END_GROUP\n')
				
				
		else:
			accepted_match = matches[int(answer)-1].split()
			jstat_resids,params,bses,chi_resids = mkl.calculate_resids([accepted_match])
			source = get_srcg(accepted_match)
			source.SI = params[0][0]
			source.intercept = params[0][1]
			source.SI_err = bses[0][0]
			source.intercept_err = bses[0][1]
			source.chi_resid = chi_resids[0]
			source.epsilon_red = jstat_resids[0]
			source.inspected = 2
			
			source.accept = accept_type
			if chi_resids[0]<=2:
				source.low_resids = 0
			else:
				source.low_resids = 1
				
			for cat,name in zip(source.cats[1:],source.names[1:]):
				if cat == 'vlssr':
					source.vlssr = name
				elif cat == 'mrc':
					source.mrc = name
				elif cat == 'sumss':
					source.sumss = name
				elif cat == 'nvss':
					source.nvss = name
				elif cat == 'atg20':
					source.atg20 = name
				else:
					pass
				
			sources.append(source)
			#g_stats.accept_type = 'spectral'
			sources_stats.append(g_stats)

for name in sorted(reject_list): reject_list_file.write(name+'\n')
#for name in sorted(v5_list): v5_list_file.write(name+'\n')
		
t=Table(masked=True)

out_t = make_the_table(t,sources,sources_stats)

out_t.write("%s_choices.vot" %options.tag,format='votable',overwrite=True)
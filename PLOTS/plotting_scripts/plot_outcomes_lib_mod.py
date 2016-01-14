#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import patches
from matplotlib.patches import Ellipse
from itertools import combinations
from astropy.wcs import WCS
from wcsaxes import WCSAxes
from sys import path
#path.append('/home/jline/Documents/cataloguing/FHD/carroll/paper/complete_v11/investigate_possibles_fix')
import make_table_lib_mod as mkl
import subprocess
#import make_modifications as mdd

dr = np.pi/180.0

##Variables that get used by most of the functions here
global closeness
global high_prob
global low_prob
global chi_thresh
global jstat_thresh
global num_freqs
global split
global matched_cats

##Set up a load of markers, colours and alpha values
markers = ['o','*','s','^','D','8','H','>','<','8','v','d']
marker_sizes = np.array([8,11,8,8,7,10,11,11,11,11,11,11,11]) + 3
marker_colours = ['m', "b", "#e5e500", "r", "g", 'k', "#c0531f",'#660066','#000099','y','#990000','#003300']
ell_colours1 = ["#C370C8", "#698ACD", "#AA6600", "#D84C77", "#5FC34F", 'k', "#D46733",'m','b','y','r','g']
ell_colours2 = ['m', "b", "#8B5A00", "r", "g", 'k', "#c0531f",'#660066','#000099','y','#990000','#003300']
alphas = [0.4,0.4,0.4,0.4,0.4,0.4,0.3,0.45,0.5,0.4,0.5,0.4]



#marker_sizes = [15,15,15,15,15,15,15,15,15,15,15,15]


##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
def plot_errors(style,colour,freq,flux,ferr,name,size,ax):
	'''Plots errorbars and markers with no line'''
	ax.errorbar(freq,flux,ferr,
	marker=style,ms=size,mfc=colour,mec=colour,ecolor=colour,markeredgewidth=1,label=name,linestyle='None')
	
def plot_errors_comb(style,colour,freq,flux,ferr,name,size,ax):
	'''Plots errorbars and markers with no line'''
	ax.errorbar(freq,flux,ferr,
	marker=style,ms=size,mfc=colour,mec='k',ecolor=colour,markeredgewidth=1,label=name,linestyle='None')
	
def plot_pos(style,colour,ra,dec,rerr,derr,name,size,ax,proj):
	'''Plots a single point with x and y erros bars'''
	if proj==1.0:
		p = ax.errorbar(ra,dec,derr,rerr,marker=style,ms=size,mfc=colour,mec=colour,ecolor=colour,markeredgewidth=1,label=name,linestyle='None')
		#p = ax.errorbar(ra,dec,derr,rerr,marker=style,ms=size,mfc='none',mec=colour,ecolor=colour,markeredgewidth=1.5,elinewidth=1,label=name,linestyle='None')
	else:
		p = ax.errorbar(ra,dec,derr,rerr,marker=style,ms=size,mfc=colour,mec=colour,ecolor=colour,markeredgewidth=1,label=name,linestyle='None',transform=proj)
		#p = ax.errorbar(ra,dec,derr,rerr,marker=style,ms=size,mfc='none',mec=colour,ecolor=colour,markeredgewidth=1.5,elinewidth=1,label=name,linestyle='None',transform=proj)
	return p

def plot_pos_comb(style,colour,ra,dec,rerr,derr,name,size,ax,proj):
	'''Plots a single point with x and y erros bars'''
	if proj==1.0:
		p = ax.errorbar(ra,dec,derr,rerr,marker=style,ms=size,mfc=colour,mec='k',ecolor=colour,markeredgewidth=1.2,label=name,linestyle='None')
	else:
		p = ax.errorbar(ra,dec,derr,rerr,marker=style,ms=size,mfc=colour,mec='k',ecolor=colour,markeredgewidth=1.2,label=name,linestyle='None',transform=proj)
	return p
	
def plt_ell(ra,dec,height,width,PA,ax,colour,colour2,alpha,proj):
	'''Plots an ellipse - either plots on the ax_main which uses a wcs projection
	or on the smaller subplots which don't need transforming'''
	##Position Angle measures angle from direction to NCP towards increasing RA (east)
	##Matplotlib plots the angle from the increasing y-axis toward DECREASING x-axis
	##so have to put in the PA as negative
	if proj==1.0:
		ell = Ellipse([ra,dec],width=width,height=height,angle=-PA,linewidth=1.5)  ##minus???
		ell.set_facecolor('none')
		ell.set_edgecolor(colour2)
		ax.add_artist(ell)
		ell = Ellipse([ra,dec],width=width,height=height,angle=-PA,linewidth=1.5)
		ell.set_facecolor(colour)
		ell.set_alpha(0.3)
		ax.add_artist(ell)
	else:
		ell = Ellipse([ra,dec],width=width,height=height,angle=-PA,transform=proj,linewidth=1.5)  ##minus???
		ell.set_facecolor('none')
		ell.set_edgecolor(colour2)
		ax.add_artist(ell)
		ell = Ellipse([ra,dec],width=width,height=height,angle=-PA,transform=proj,linewidth=1.5)
		ell.set_facecolor(colour)
		ell.set_alpha(0.3)
		ax.add_artist(ell)
	
##POSSIBLE EXTENSION - MAKE GENERIC SO IT CYCLES THROUGH SOME COLOURS, NOT SPECIFIED COLOURS
##FOR A PARTICULAR CATALOGUE
def plot_all(cat,name,ra,rerr,dec,derr,major,minor,PA,ax,proj):
	##Plot the colour by index of catalogue in matched_cats 
	ind = matched_cats.index(cat)
	plot_pos(markers[ind],marker_colours[ind],ra,dec,rerr,derr,name,marker_sizes[ind],ax,proj)
	if cat=='mrc':
		pass
	else:
		if float(minor)!=-100000.0:
			if float(major)!=-100000.0:
				plt_ell(ra,dec,float(major),float(minor),float(PA),ax,ell_colours1[ind],ell_colours2[ind],alphas[ind],proj)
				
##--------------------------------------------------------------------------------------------------------------------

def get_allinfo(all_info,good_names):
	'''Takes a list of strings. Each string is a line containing the information for a single source
	in a matched group from an output file of calculate_bayes.py. Gets all of the information
	from each entry and returns them in a source_group() class'''
	src_all = mkl.source_group()
	for entry in all_info:
		info = entry.split()
		if info[1] in good_names or all_info.index(entry) == 0 or good_names==['all']:
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
				src_all.freqs.append(np.array([float(info[6])]))          
				src_all.fluxs.append(np.array([float(info[7])]))
				src_all.ferrs.append(np.array([float(info[8])]))
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
					if np.isnan(float(info[7+(3*i)])) == False or float(info[7+(3*i)]) != -100000.0:
						freqs.append(float(info[6+(3*i)]))
						fluxs.append(float(info[7+(3*i)]))
						ferrs.append(float(info[8+(3*i)]))
					else:
						freqs.append(-100000.0)
						fluxs.append(-100000.0)
						ferrs.append(-100000.0)
				src_all.freqs.append(np.array(freqs))
				src_all.fluxs.append(np.array(fluxs))
				src_all.ferrs.append(np.array(ferrs))
					#src_all.freqs.append(float(info[6+(3*i)]))
					#src_all.fluxs.append(float(info[7+(3*i)]))
					#src_all.ferrs.append(float(info[8+(3*i)]))
	return src_all

def get_srcg(src_all,info,good_names):
	'''Takes a string which contains the information for all sources in particular combination
	Uses num_freqs to work out where each piece of information is, and then return the relevant
	information in a source_group() class. Any catalogues with no matches are entered as -10000.0'''
	src_g = mkl.source_group()
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
		
		#print info[ind+1],good_names
		
		cat_slot = matched_cats[j]
		
		if cat_slot in src_all.cats:
			srcall_ind = src_all.cats.index(cat_slot)
			freqss = []
			fluxss = []
			ferrss = []
			for k in xrange(num_freq):
				##Test to see if flux is a nan or -100000.0; make sure all flux/freq info is -100000.0 if so
				if np.isnan(src_all.fluxs[srcall_ind][k]) == False or float(src_all.fluxs[srcall_ind][k]) != -100000.0:
					freqss.append(float(src_all.freqs[srcall_ind][k]))
					fluxss.append(float(src_all.fluxs[srcall_ind][k]))
					ferrss.append(float(src_all.ferrs[srcall_ind][k]))
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
			src_g.cats.append(src_all.cats[srcall_ind])
			src_g.names.append(src_all.names[srcall_ind])
			src_g.ras.append(src_all.ras[srcall_ind])
			src_g.rerrs.append(src_all.rerrs[srcall_ind])
			src_g.decs.append(src_all.decs[srcall_ind])
			src_g.derrs.append(src_all.derrs[srcall_ind])
			src_g.majors.append(src_all.majors[srcall_ind])
			src_g.minors.append(src_all.minors[srcall_ind])
			src_g.PAs.append(src_all.PAs[srcall_ind])
			src_g.flags.append(src_all.flags[srcall_ind])
			src_g.IDs.append(src_all.IDs[srcall_ind])
			src_g.prob = float(info[-1])
			#freqss = []
			#fluxss = []
			#ferrss = []
			#for k in xrange(num_freq):
				###Test to see if flux is a nan or -100000.0; make sure all flux/freq info is -100000.0 if so
				#if np.isnan(float(info[7+ind+(3*k)])) == False or float(info[7+(3*i)]) != -100000.0:
					#freqss.append(float(info[6+ind+(3*k)]))
					#fluxss.append(float(info[7+ind+(3*k)]))
					#ferrss.append(float(info[8+ind+(3*k)]))
				#else:
					#freqss.append(-100000.0)
					#fluxss.append(-100000.0)
					#ferrss.append(-100000.0)
				##if float(info[6+ind+(3*k)])!=-100000.0: all_freqs.append(float(info[6+ind+(3*k)]))
				##if float(info[7+ind+(3*k)])!=-100000.0: all_fluxs.append(float(info[7+ind+(3*k)]))
				##if float(info[8+ind+(3*k)])!=-100000.0: all_ferrs.append(float(info[8+ind+(3*k)]))
			#src_g.freqs.append(freqss)
			#src_g.fluxs.append(fluxss)
			#src_g.ferrs.append(ferrss)
			#src_g.cats.append(info[ind])
			#src_g.names.append(info[ind+1])
			#src_g.ras.append(float(info[ind+2]))
			#src_g.rerrs.append(float(info[ind+3]))
			#src_g.decs.append(float(info[ind+4]))
			#src_g.derrs.append(float(info[ind+5]))
			#src_g.majors.append(float(info[ind+9+((num_freq-1)*3)]))
			#src_g.minors.append(float(info[ind+10+((num_freq-1)*3)]))
			#src_g.PAs.append(float(info[ind+11+((num_freq-1)*3)]))
			#src_g.flags.append(info[ind+12+((num_freq-1)*3)])
			#src_g.IDs.append(info[ind+13+((num_freq-1)*3)])
			#src_g.prob = float(info[-1])
		else:
			freqss = []
			fluxss = []
			ferrss = []
			for k in xrange(num_freq):
				freqss.append(-100000.0)
				fluxss.append(-100000.0)
				ferrss.append(-100000.0)
				#if float(info[6+ind+(3*k)])!=-100000.0: all_freqs.append(float(info[6+ind+(3*k)]))
				#if float(info[7+ind+(3*k)])!=-100000.0: all_fluxs.append(float(info[7+ind+(3*k)]))
				#if float(info[8+ind+(3*k)])!=-100000.0: all_ferrs.append(float(info[8+ind+(3*k)]))
			src_g.freqs.append(freqss)
			src_g.fluxs.append(fluxss)
			src_g.ferrs.append(ferrss)
			src_g.cats.append('-100000.0')
			src_g.names.append('-100000.0')
			src_g.ras.append(-100000.0)
			src_g.rerrs.append(-100000.0)
			src_g.decs.append(-100000.0)
			src_g.derrs.append(-100000.0)
			src_g.majors.append(-100000.0)
			src_g.minors.append(-100000.0)
			src_g.PAs.append(-100000.0)
			src_g.flags.append('-100000.0')
			src_g.IDs.append('-100000.0')
			src_g.prob = float(info[-1])
	return src_g


#def plot_ind(match,ax,ind_ax,ax_spectral,ra_bottom,ra_top,dec_bottom,dec_top,dom_crit,comb_crit):
def plot_ind(match,ind_ax,ax_spectral,dom_crit,comb_crit):
	'''Takes a string of information of a particular combination and uses it
	to create a plot of the single combination, and fit and plot a line to
	the spectral information. Returns the positional probability, fitted paramaters
	and the residuals of the fit'''
	
	##Get the information from the particular match given
	info = match.split()
	indexes = [(14+((i-1)*3)) for i in num_freqs]
	starts = [0]
	for i in xrange(len(indexes)-1): starts.append(sum(indexes[:i+1]))
	fluxs = []
	freqs = []
	ferrs = []
	for j in xrange(len(starts)): 
		ind = starts[j]
		cat = info[ind]
		if cat!='-100000.0':
			#num_freq = num_freqs[j]
			freq = num_freqs[j]
			name = info[ind+1]
			ra = float(info[ind+2])
			rerr = float(info[ind+3])
			dec = float(info[ind+4])
			derr = float(info[ind+5])
			nu = float(info[ind+6])
			#flux = float(info[ind+7])
			#ferr = float(info[ind+8])/flux
			for k in xrange(freq):
				if info[7+ind+(3*k)]!='-100000.0':
					if np.isnan(float(info[7+ind+(3*k)])) == False:
						freqs.append(float(info[6+ind+(3*k)]))
						fluxs.append(float(info[7+ind+(3*k)]))
						ferrs.append(float(info[8+ind+(3*k)])/float(info[7+ind+(3*k)]))
			major = info[ind+9+((freq-1)*3)]
			minor = info[ind+10+((freq-1)*3)]
			PA = info[ind+11+((freq-1)*3)]
			#fluxs.append(flux)
			#freqs.append(nu)
			#ferrs.append(ferr)
			##Plot each source on the individual combo plot
			#plot_all(cat,name,ra,rerr,dec,derr,major,minor,PA,ax,1.0)
	##Sort the frequencies, fluxes and log them
	log_fluxs = np.log([flux for (freq,flux) in sorted(zip(freqs,fluxs),key=lambda pair: pair[0])])
	sorted_ferrs = np.array([ferr for (freq,ferr) in sorted(zip(freqs,ferrs),key=lambda pair: pair[0])])
	log_freqs = np.log(sorted(freqs))
	ferrs = np.array(ferrs)
	prob = info[-1]
	
	##Fit a line using weighted least squares and plot it
	lin_fit,jstat,bse,chi_red = mkl.fit_line(log_freqs,log_fluxs,sorted_ferrs)
	
	#ax.text(0.5,0.925, 'P$_{%d}$=%.2f' %(ind_ax+1,float(prob)), transform=ax.transAxes,verticalalignment='center',horizontalalignment='center',fontsize=16)
	#ax.text(0.5,0.06, '$\epsilon_{%d}$=%.2f $\chi_{%d}$=%.1f' %(ind_ax+1,jstat,ind_ax+1,chi_red),
		#transform=ax.transAxes,verticalalignment='center',horizontalalignment='center',fontsize=16)
	
	#ax.set_xlim(ra_bottom,ra_top)
	#ax.set_ylim(dec_bottom,dec_top)
	
	###Plot RA backwards
	#ax.invert_xaxis()
	
	##Plot the fitted spectral line, and return the plot object so we can create a legend from it
	if dom_crit=='No dom. source':
		if 'split' in comb_crit:
			spec_plot = 'na'
		else:
			spec_plot, = ax_spectral.plot(np.exp(log_freqs),np.exp(lin_fit.fittedvalues),linestyle='-',linewidth=1,alpha=0.3)
	else:
		spec_plot, = ax_spectral.plot(np.exp(log_freqs),np.exp(lin_fit.fittedvalues),linestyle='-',linewidth=1,alpha=0.7)
	
	return prob,jstat,spec_plot,lin_fit.params

def make_left_plots(fig,main_dims,spec_dims,ra_main,dec_main):
	
	##A fits image header with which to create a wcs with
	header = { 'NAXIS'  : 2,             ##Number of data axis
    'NAXIS1' : 10,                  ##Length of X axis
    'CTYPE1' : 'RA---SIN',           ##Projection type of X axis
	'CRVAL1' : ra_main,        ##Central X world coord value
	'CRPIX1' : 5,                    ##Central X Pixel value
	'CUNIT1' : 'deg',                ##Unit of X axes
	'CDELT1' : -0.001*np.cos(dec_main*(np.pi/180.0)),              ##Size of pixel in world co-ord
	'NAXIS2' : 10,                  ##Length of X axis
	'CTYPE2' : 'DEC--SIN',           ##Projection along Y axis
	'CRVAL2' : dec_main,                   ##Central Y world coord value
	'CRPIX2' : 5,                    ##Central Y Pixel value
	'CUNIT2' : 'deg',                ##Unit of Y world coord
	'CDELT2' : +0.001      		     ##Size of pixel in deg
	} 
	
	##Create the ws, and the main axis based on that. Plot top left
	wcs = WCS(header=header)
	ax_main = WCSAxes(fig, main_dims, wcs=wcs)
	fig.add_axes(ax_main)
	tr_fk5 = ax_main.get_transform("fk5")
	
	#ax_main.set_title("All sources within 3'.0")
	ax_main.text(0.01,0.93,"All sources within search area",verticalalignment='bottom',horizontalalignment='left', transform=ax_main.transAxes,fontsize=16)
	
	##Create bottom left plot with log-log axes - set the error bars to plot
	##even if they go off the edge of the plot
	ax_spectral = fig.add_axes(spec_dims)
	ax_spectral.set_xscale("log",nonposx='clip')
	ax_spectral.set_yscale("log",nonposy='clip')
	
	return ax_main,ax_spectral,tr_fk5,wcs

def make_right_plots(fig,main_dims,spec_dims,ra_main,dec_main):
	
	##A fits image header with which to create a wcs with
	header = { 'NAXIS'  : 2,             ##Number of data axis
    'NAXIS1' : 10,                  ##Length of X axis
    'CTYPE1' : 'RA---SIN',           ##Projection type of X axis
	'CRVAL1' : ra_main,        ##Central X world coord value
	'CRPIX1' : 5,                    ##Central X Pixel value
	'CUNIT1' : 'deg',                ##Unit of X axes
	'CDELT1' : -0.001*np.cos(dec_main*(np.pi/180.0)),              ##Size of pixel in world co-ord
	'NAXIS2' : 10,                  ##Length of X axis
	'CTYPE2' : 'DEC--SIN',           ##Projection along Y axis
	'CRVAL2' : dec_main,                   ##Central Y world coord value
	'CRPIX2' : 5,                    ##Central Y Pixel value
	'CUNIT2' : 'deg',                ##Unit of Y world coord
	'CDELT2' : +0.001      		     ##Size of pixel in deg
	} 
	
	##Create the ws, and the main axis based on that. Plot top left
	wcs = WCS(header=header)
	ax_main = WCSAxes(fig, main_dims, wcs=wcs)
	fig.add_axes(ax_main)
	tr_fk5 = ax_main.get_transform("fk5")
	
	#ax_main.set_title("All sources within 3'.0")
	ax_main.text(0.01,0.93,"Sources included by user",verticalalignment='bottom',horizontalalignment='left', transform=ax_main.transAxes,fontsize=16)
	
	##Create bottom left plot with log-log axes - set the error bars to plot
	##even if they go off the edge of the plot
	ax_spectral = fig.add_axes(spec_dims)
	ax_spectral.set_xscale("log",nonposx='clip')
	ax_spectral.set_yscale("log",nonposy='clip')
	
	return ax_main,ax_spectral,tr_fk5,wcs

def fill_left_plots(all_info,ra_main,dec_main,ax_main,ax_spectral,tr_fk5,wcs,all_fluxs,ra_down_lim,ra_up_lim,dec_down_lim,dec_up_lim,delta_RA):
	'''Get the information and plot the positions and fluxs on the left hand side plots'''
	ras = []
	decs = []
	all_freqs = []
	all_ferrs = []
	for i in xrange(len(all_info)):
		info=all_info[i].split()
		cat = info[0]
		name = info[1]
		ra = float(info[2])
		ras.append(ra)
		rerr = float(info[3])
		dec = float(info[4])
		decs.append(dec)
		derr = float(info[5])
		major = info[-5]
		minor = info[-4]
		PA = info[-3]
		ID = info[-1]
		##If base catalogue, plot error and mathcing ellipses
		if i==0:
			##Calculate the change in RA for half the user given resolution (and constant dec) using spherical trigo
			error_RA = np.arccos((np.cos(closeness*dr)-np.sin(dec*dr)**2)/np.cos(dec*dr)**2)/dr
			
			##Plot an error ellipse of the base cat error + resolution
			ell = patches.Ellipse((ra_main,dec_main),2*(rerr+error_RA),2*(derr+closeness),angle=0,
				transform=tr_fk5,linestyle='dashed',fc='none',lw=1.1,color='gray')
			ax_main.add_patch(ell)
			##Plot a circle of the match radius
			ell = patches.Ellipse((ra_main,dec_main),2*delta_RA,4*(closeness),angle=0,
				transform=tr_fk5,linestyle='dashdot',fc='none',lw=1.1,color='k')
			ax_main.add_patch(ell)
		
		##Plot positions and elliptical fits
		plot_all(cat,name,ra,rerr,dec,derr,major,minor,PA,ax_main,tr_fk5)
		
		##See if one or more flux for a source, and plot fluxes with errorbars
		cat_ind = matched_cats.index(cat)
		if len(info)==14:
			freq = float(info[6])
			flux = float(info[7])
			ferr = float(info[8])
			plot_errors(markers[cat_ind],marker_colours[cat_ind],freq,flux,ferr,name,marker_sizes[cat_ind],ax_spectral)
			all_fluxs.append(flux)
			all_freqs.append(freq)
			all_ferrs.append(ferr)
		##If the catalogue has more than one frequency:
		else:
			extra = (len(info)-14) / 3
			freqs = []
			fluxs = []
			ferrs = []
			for i in xrange(extra+1):
				if info[7+(3*i)]!='-100000.0' and info[7+(3*i)]!='nan':
					freqs.append(info[6+(3*i)])
					fluxs.append(info[7+(3*i)])
					ferrs.append(info[8+(3*i)])
			
			for i in xrange(len(fluxs)):
				all_fluxs.append(float(fluxs[i]))
				plot_errors(markers[cat_ind],marker_colours[cat_ind],float(freqs[i]),float(fluxs[i]),float(ferrs[i]),'%s-%.1fMHz' %(name,float(freqs[i])),marker_sizes[cat_ind],ax_spectral)
			for freq in freqs: all_freqs.append(float(freq))
			for ferr in ferrs: all_ferrs.append(float(ferr))
				
	##Add some labels and coord formatting to ax_main
	ra_ax = ax_main.coords[0]
	dec_ax = ax_main.coords[1]
	ra_ax.set_axislabel('RAJ2000')
	dec_ax.set_axislabel('DECJ2000')
	ra_ax.set_major_formatter('hh:mm:ss')
	dec_ax.set_major_formatter('dd:mm:ss')
	
	##Convert axes limits to ax_main wcs, and apply
	ra_low = wcs.wcs_world2pix(ra_down_lim,dec_main,0)  ##The zero is for the orgin point of the image
	ra_high = wcs.wcs_world2pix(ra_up_lim,dec_main,0)
	dec_low = wcs.wcs_world2pix(ra_main,dec_down_lim,0)
	dec_high = wcs.wcs_world2pix(ra_main,dec_up_lim,0)
	ax_main.set_ylim(dec_low[1],dec_high[1])
	ax_main.set_xlim(ra_high[0],ra_low[0])

	##Make the labels on ax_spectral print in MHz and Jy
	max_lflux = np.log(max(all_fluxs))
	min_lflux = np.log(min(all_fluxs))
	
	freq_ticks = [freq for freq in sorted(set(all_freqs))]
	flux_ticks = np.exp(np.arange(min_lflux,max_lflux+abs(max_lflux-min_lflux)/5,abs(max_lflux-min_lflux)/5))
	ax_spectral.set_xticks(freq_ticks,minor=False)
	ax_spectral.set_xticklabels(freq_ticks)
	ax_spectral.set_yticks(flux_ticks,minor=False)
	ax_spectral.set_yticklabels(['%.3f' %flux for flux in list(flux_ticks)],fontsize=14.0)

	##Set some limits on the spextral axis - need to do it in logspace,
	##and set the edge gaps based on the data
	freq_delta = abs(np.log10(max(all_freqs))-np.log10(min(all_freqs)))/10.0
	flux_delta = abs(np.log10(max(all_fluxs))-np.log10(min(all_fluxs)))/10.0
	
	freq_max = 10**(np.log10(max(all_freqs))+freq_delta)
	freq_min = 10**(np.log10(min(all_freqs))-freq_delta)
	flux_max = 10**(np.log10(max(all_fluxs))+flux_delta)
	flux_min = 10**(np.log10(min(all_fluxs))-flux_delta)
	
	ax_spectral.set_xlim(freq_min,freq_max)
	ax_spectral.set_ylim(flux_min,flux_max)

	##Stick some grid stuff on the log log plot
	ten_steps = []
	for arr in [np.array([1e-5,2e-5,3e-5,4e-5,5e-5,6e-5,7e-5,8e-5,9e-5])*10**x for x in xrange(10)]:
		for i in list(arr): ten_steps.append(i)

	ax_spectral.set_xticks([step for step in ten_steps if step>freq_min and step<freq_max],minor=True)
	ax_spectral.set_yticks([step for step in ten_steps if step>flux_min and step<flux_max],minor=True)

	ax_spectral.xaxis.grid(True, which='minor',linestyle='dashed',alpha=0.1)
	ax_spectral.yaxis.grid(True, which='minor',linestyle='dashed',alpha=0.1)
	ax_spectral.set_xlabel(r'log$_{10}$(Frequency) (MHz)',fontsize=14)
	ax_spectral.set_ylabel(r'log$_{10}$(Flux) (Jy)',fontsize=14)
	
def fill_right_plots(all_info,ra_main,dec_main,ax_main,ax_spectral,tr_fk5,wcs,all_fluxs,ra_down_lim,ra_up_lim,dec_down_lim,dec_up_lim,delta_RA,good_names):
	'''Get the information and plot the positions and fluxs on the left hand side plots'''
	ras = []
	decs = []
	all_freqs = []
	all_ferrs = []
	for i in xrange(len(all_info)):
		info=all_info[i].split()
		cat = info[0]
		name = info[1]
		if name in good_names or good_names==['all'] or matched_cats.index(cat) == 0:
			
			
			ra = float(info[2])
			ras.append(ra)
			rerr = float(info[3])
			dec = float(info[4])
			decs.append(dec)
			derr = float(info[5])
			major = info[-5]
			minor = info[-4]
			PA = info[-3]
			ID = info[-1]
			##If base catalogue, plot error and mathcing ellipses
			if i==0:
				##Calculate the change in RA for half the user given resolution (and constant dec) using spherical trigo
				error_RA = np.arccos((np.cos(closeness*dr)-np.sin(dec*dr)**2)/np.cos(dec*dr)**2)/dr
				
				##Plot an error ellipse of the base cat error + resolution
				ell = patches.Ellipse((ra_main,dec_main),2*(rerr+error_RA),2*(derr+closeness),angle=0,
					transform=tr_fk5,linestyle='dashed',fc='none',lw=1.1,color='gray')
				ax_main.add_patch(ell)
				##Plot a circle of the match radius
				ell = patches.Ellipse((ra_main,dec_main),2*delta_RA,4*(closeness),angle=0,
					transform=tr_fk5,linestyle='dashdot',fc='none',lw=1.1,color='k')
				ax_main.add_patch(ell)
			
			##Plot positions and elliptical fits
			plot_all(cat,name,ra,rerr,dec,derr,major,minor,PA,ax_main,tr_fk5)
			
			##See if one or more flux for a source, and plot fluxes with errorbars
			cat_ind = matched_cats.index(cat)
			if len(info)==14:
				freq = float(info[6])
				flux = float(info[7])
				ferr = float(info[8])
				plot_errors(markers[cat_ind],marker_colours[cat_ind],freq,flux,ferr,name,marker_sizes[cat_ind],ax_spectral)
				all_fluxs.append(flux)
				all_freqs.append(freq)
				all_ferrs.append(ferr)
			##If the catalogue has more than one frequency:
			else:
				extra = (len(info)-14) / 3
				freqs = []
				fluxs = []
				ferrs = []
				for i in xrange(extra+1):
					if info[7+(3*i)]!='-100000.0' and info[7+(3*i)]!='nan':
						freqs.append(info[6+(3*i)])
						fluxs.append(info[7+(3*i)])
						ferrs.append(info[8+(3*i)])
				
				for i in xrange(len(fluxs)):
					all_fluxs.append(float(fluxs[i]))
					plot_errors(markers[cat_ind],marker_colours[cat_ind],float(freqs[i]),float(fluxs[i]),float(ferrs[i]),'%s-%.1fMHz' %(name,float(freqs[i])),marker_sizes[cat_ind],ax_spectral)
				for freq in freqs: all_freqs.append(float(freq))
				for ferr in ferrs: all_ferrs.append(float(ferr))
	
	##Add some labels and coord formatting to ax_main
	ra_ax = ax_main.coords[0]
	dec_ax = ax_main.coords[1]
	ra_ax.set_axislabel('RAJ2000')
	dec_ax.set_axislabel('DECJ2000')
	ra_ax.set_major_formatter('hh:mm:ss')
	dec_ax.set_major_formatter('dd:mm:ss')
	
	dec_ax.set_ticklabel_position('r')
	dec_ax.set_axislabel_position('r')
	
	
	
	##Convert axes limits to ax_main wcs, and apply
	ra_low = wcs.wcs_world2pix(ra_down_lim,dec_main,0)  ##The zero is for the orgin point of the image
	ra_high = wcs.wcs_world2pix(ra_up_lim,dec_main,0)
	dec_low = wcs.wcs_world2pix(ra_main,dec_down_lim,0)
	dec_high = wcs.wcs_world2pix(ra_main,dec_up_lim,0)
	ax_main.set_ylim(dec_low[1],dec_high[1])
	ax_main.set_xlim(ra_high[0],ra_low[0])

	##Make the labels on ax_spectral print in MHz and Jy
	max_lflux = np.log(max(all_fluxs))
	min_lflux = np.log(min(all_fluxs))
	
	freq_ticks = [freq for freq in sorted(set(all_freqs))]
	flux_ticks = np.exp(np.arange(min_lflux,max_lflux+abs(max_lflux-min_lflux)/5,abs(max_lflux-min_lflux)/5))
	ax_spectral.set_xticks(freq_ticks,minor=False)
	ax_spectral.set_xticklabels(freq_ticks)
	ax_spectral.set_yticks(flux_ticks,minor=False)
	ax_spectral.set_yticklabels(['%.3f' %flux for flux in list(flux_ticks)],fontsize=14.0)

	##Set some limits on the spextral axis - need to do it in logspace,
	##and set the edge gaps based on the data
	freq_delta = abs(np.log10(max(all_freqs))-np.log10(min(all_freqs)))/10.0
	flux_delta = abs(np.log10(max(all_fluxs))-np.log10(min(all_fluxs)))/10.0
	
	freq_max = 10**(np.log10(max(all_freqs))+freq_delta)
	freq_min = 10**(np.log10(min(all_freqs))-freq_delta)
	flux_max = 10**(np.log10(max(all_fluxs))+flux_delta)
	flux_min = 10**(np.log10(min(all_fluxs))-flux_delta)
	
	ax_spectral.set_xlim(freq_min,freq_max)
	ax_spectral.set_ylim(flux_min,flux_max)

	##Stick some grid stuff on the log log plot
	ten_steps = []
	for arr in [np.array([1e-5,2e-5,3e-5,4e-5,5e-5,6e-5,7e-5,8e-5,9e-5])*10**x for x in xrange(10)]:
		for i in list(arr): ten_steps.append(i)
		
	ax_spectral.yaxis.tick_right()
	ax_spectral.yaxis.set_label_position("right")

	ax_spectral.set_xticks([step for step in ten_steps if step>freq_min and step<freq_max],minor=True)
	ax_spectral.set_yticks([step for step in ten_steps if step>flux_min and step<flux_max],minor=True)

	ax_spectral.xaxis.grid(True, which='minor',linestyle='dashed',alpha=0.1)
	ax_spectral.yaxis.grid(True, which='minor',linestyle='dashed',alpha=0.1)
	ax_spectral.set_xlabel(r'log$_{10}$(Frequency) (MHz)',fontsize=14)
	
	
	
	ax_spectral.set_ylabel(r'log$_{10}$(Flux) (Jy)',fontsize=14)
	
def fill_right_plots_nosource(all_info,ra_main,dec_main,ax_main,ax_spectral,tr_fk5,wcs,all_fluxs,ra_down_lim,ra_up_lim,dec_down_lim,dec_up_lim,delta_RA):
	'''Get the information and plot the positions and fluxs on the left hand side plots'''
	ras = []
	decs = []
	all_freqs = []
	all_ferrs = []
	for i in xrange(len(all_info)):
		info=all_info[i].split()
		cat = info[0]
		name = info[1]
		ra = float(info[2])
		ras.append(ra)
		rerr = float(info[3])
		dec = float(info[4])
		decs.append(dec)
		derr = float(info[5])
		major = info[-5]
		minor = info[-4]
		PA = info[-3]
		ID = info[-1]
		##If base catalogue, plot error and mathcing ellipses
		if i==0:
			##Calculate the change in RA for half the user given resolution (and constant dec) using spherical trigo
			error_RA = np.arccos((np.cos(closeness*dr)-np.sin(dec*dr)**2)/np.cos(dec*dr)**2)/dr
			
			##Plot an error ellipse of the base cat error + resolution
			ell = patches.Ellipse((ra_main,dec_main),2*(rerr+error_RA),2*(derr+closeness),angle=0,
				transform=tr_fk5,linestyle='dashed',fc='none',lw=1.1,color='gray')
			ax_main.add_patch(ell)
			##Plot a circle of the match radius
			ell = patches.Ellipse((ra_main,dec_main),2*delta_RA,4*(closeness),angle=0,
				transform=tr_fk5,linestyle='dashdot',fc='none',lw=1.1,color='k')
			ax_main.add_patch(ell)
		
		##Plot positions and elliptical fits
		plot_all(cat,name,ra,rerr,dec,derr,major,minor,PA,ax_main,tr_fk5)
		
		##See if one or more flux for a source, and plot fluxes with errorbars
		cat_ind = matched_cats.index(cat)
		if len(info)==14:
			freq = float(info[6])
			flux = float(info[7])
			ferr = float(info[8])
			if cat_ind!=0: plot_errors(markers[cat_ind],marker_colours[cat_ind],freq,flux,ferr,name,marker_sizes[cat_ind],ax_spectral)
			all_fluxs.append(flux)
			all_freqs.append(freq)
			all_ferrs.append(ferr)
		##If the catalogue has more than one frequency:
		else:
			extra = (len(info)-14) / 3
			freqs = []
			fluxs = []
			ferrs = []
			for i in xrange(extra+1):
				freqs.append(info[6+(3*i)])
				fluxs.append(info[7+(3*i)])
				ferrs.append(info[8+(3*i)])
			for flux in fluxs: 
				all_fluxs.append(float(flux))
				if float(flux)!=100000.0 and cat_ind!=0: plot_errors(markers[cat_ind],marker_colours[cat_ind],float(freq),float(flux),float(ferr),'%s-%.1fMHz' %(name,freq),marker_sizes[cat_ind],ax_spectral)
			for freq in freqs: all_freqs.append(float(freq))
			for ferr in ferrs: all_ferrs.append(float(ferr))
	
	##Add some labels and coord formatting to ax_main
	ra_ax = ax_main.coords[0]
	dec_ax = ax_main.coords[1]
	ra_ax.set_axislabel('RAJ2000')
	dec_ax.set_axislabel('DECJ2000')
	ra_ax.set_major_formatter('hh:mm:ss')
	dec_ax.set_major_formatter('dd:mm:ss')
	
	dec_ax.set_ticklabel_position('r')
	dec_ax.set_axislabel_position('r')
	
	
	
	##Convert axes limits to ax_main wcs, and apply
	ra_low = wcs.wcs_world2pix(ra_down_lim,dec_main,0)  ##The zero is for the orgin point of the image
	ra_high = wcs.wcs_world2pix(ra_up_lim,dec_main,0)
	dec_low = wcs.wcs_world2pix(ra_main,dec_down_lim,0)
	dec_high = wcs.wcs_world2pix(ra_main,dec_up_lim,0)
	ax_main.set_ylim(dec_low[1],dec_high[1])
	ax_main.set_xlim(ra_high[0],ra_low[0])

	##Make the labels on ax_spectral print in MHz and Jy
	max_lflux = np.log(max(all_fluxs))
	min_lflux = np.log(min(all_fluxs))
	
	freq_ticks = [freq for freq in sorted(set(all_freqs))]
	flux_ticks = np.exp(np.arange(min_lflux,max_lflux+abs(max_lflux-min_lflux)/5,abs(max_lflux-min_lflux)/5))
	ax_spectral.set_xticks(freq_ticks,minor=False)
	ax_spectral.set_xticklabels(freq_ticks)
	ax_spectral.set_yticks(flux_ticks,minor=False)
	ax_spectral.set_yticklabels(['%.3f' %flux for flux in list(flux_ticks)],fontsize=14.0)

	##Set some limits on the spextral axis - need to do it in logspace,
	##and set the edge gaps based on the data
	freq_delta = abs(np.log10(max(all_freqs))-np.log10(min(all_freqs)))/10.0
	flux_delta = abs(np.log10(max(all_fluxs))-np.log10(min(all_fluxs)))/10.0
	
	freq_max = 10**(np.log10(max(all_freqs))+freq_delta)
	freq_min = 10**(np.log10(min(all_freqs))-freq_delta)
	flux_max = 10**(np.log10(max(all_fluxs))+flux_delta)
	flux_min = 10**(np.log10(min(all_fluxs))-flux_delta)
	
	ax_spectral.set_xlim(freq_min,freq_max)
	ax_spectral.set_ylim(flux_min,flux_max)

	##Stick some grid stuff on the log log plot
	ten_steps = []
	for arr in [np.array([1e-5,2e-5,3e-5,4e-5,5e-5,6e-5,7e-5,8e-5,9e-5])*10**x for x in xrange(10)]:
		for i in list(arr): ten_steps.append(i)
		
	ax_spectral.yaxis.tick_right()
	ax_spectral.yaxis.set_label_position("right")

	ax_spectral.set_xticks([step for step in ten_steps if step>freq_min and step<freq_max],minor=True)
	ax_spectral.set_yticks([step for step in ten_steps if step>flux_min and step<flux_max],minor=True)

	ax_spectral.xaxis.grid(True, which='minor',linestyle='dashed',alpha=0.1)
	ax_spectral.yaxis.grid(True, which='minor',linestyle='dashed',alpha=0.1)
	ax_spectral.set_xlabel(r'log$_{10}$(Frequency) (MHz)',fontsize=14)
	
	
	
	ax_spectral.set_ylabel(r'log$_{10}$(Flux) (Jy)',fontsize=14)
	
	
def fill_left_plots_reduced(all_info,ra_main,dec_main,ax_main,ax_spectral,tr_fk5,wcs,all_fluxs,ra_down_lim,ra_up_lim,dec_down_lim,dec_up_lim,delta_RA):
	'''Get the information and plot the positions and fluxs on the left hand side plots'''
	ras = []
	decs = []
	all_freqs = []
	all_ferrs = []
	for i in xrange(len(all_info)):
		info=all_info[i].split()
		cat = info[0]
		name = info[1]
		ra = float(info[2])
		ras.append(ra)
		rerr = float(info[3])
		dec = float(info[4])
		decs.append(dec)
		derr = float(info[5])
		major = info[-5]
		minor = info[-4]
		PA = info[-3]
		ID = info[-1]
		##If base catalogue, plot error and mathcing ellipses
		if i==0:
			##Calculate the change in RA for half the user given resolution (and constant dec) using spherical trigo
			error_RA = np.arccos((np.cos(closeness*dr)-np.sin(dec*dr)**2)/np.cos(dec*dr)**2)/dr
			
			##Plot an error ellipse of the base cat error + resolution
			ell = patches.Ellipse((ra_main,dec_main),2*(rerr+error_RA),2*(derr+closeness),angle=0,
				transform=tr_fk5,linestyle='dashed',fc='none',lw=1.1,color='gray')
			ax_main.add_patch(ell)
			##Plot a circle of the match radius
			ell = patches.Ellipse((ra_main,dec_main),2*delta_RA,4*(closeness),angle=0,
				transform=tr_fk5,linestyle='dashdot',fc='none',lw=1.1,color='k')
			ax_main.add_patch(ell)
		
		##Plot positions and elliptical fits
		plot_all(cat,name,ra,rerr,dec,derr,major,minor,PA,ax_main,tr_fk5)
		
		##See if one or more flux for a source, and plot fluxes with errorbars
		cat_ind = matched_cats.index(cat)
		if len(info)==14:
			freq = float(info[6])
			flux = float(info[7])
			ferr = float(info[8])
			plot_errors(markers[cat_ind],marker_colours[cat_ind],freq,flux,ferr,name,marker_sizes[cat_ind],ax_spectral)
			all_fluxs.append(flux)
			all_freqs.append(freq)
			all_ferrs.append(ferr)
		##If the catalogue has more than one frequency:
		else:
			extra = (len(info)-14) / 3
			freqs = []
			fluxs = []
			ferrs = []
			for i in xrange(extra+1):
				freqs.append(info[6+(3*i)])
				fluxs.append(info[7+(3*i)])
				ferrs.append(info[8+(3*i)])
			for flux in fluxs: 
				all_fluxs.append(float(flux))
				if float(flux)!=100000.0: plot_errors(markers[cat_ind],marker_colours[cat_ind],float(freq),float(flux),float(ferr),'%s-%.1fMHz' %(name,freq),marker_sizes[cat_ind],ax_spectral)
			for freq in freqs: all_freqs.append(float(freq))
			for ferr in ferrs: all_ferrs.append(float(ferr))
				
	##Add some labels and coord formatting to ax_main
	ra_ax = ax_main.coords[0]
	dec_ax = ax_main.coords[1]
	ra_ax.set_axislabel('RAJ2000')
	dec_ax.set_axislabel('DECJ2000')
	ra_ax.set_major_formatter('hh:mm:ss')
	dec_ax.set_major_formatter('dd:mm:ss')
	
	##Convert axes limits to ax_main wcs, and apply
	ra_low = wcs.wcs_world2pix(ra_down_lim,dec_main,0)  ##The zero is for the orgin point of the image
	ra_high = wcs.wcs_world2pix(ra_up_lim,dec_main,0)
	dec_low = wcs.wcs_world2pix(ra_main,dec_down_lim,0)
	dec_high = wcs.wcs_world2pix(ra_main,dec_up_lim,0)
	ax_main.set_ylim(dec_low[1],dec_high[1])
	ax_main.set_xlim(ra_high[0],ra_low[0])

	##Make the labels on ax_spectral print in MHz and Jy
	max_lflux = np.log(max(all_fluxs))
	min_lflux = np.log(min(all_fluxs))
	
	freq_ticks = [freq for freq in sorted(set(all_freqs))]
	flux_ticks = np.exp(np.arange(min_lflux,max_lflux+abs(max_lflux-min_lflux)/5,abs(max_lflux-min_lflux)/5))
	ax_spectral.set_xticks(freq_ticks,minor=False)
	ax_spectral.set_xticklabels(freq_ticks)
	ax_spectral.set_yticks(flux_ticks,minor=False)
	ax_spectral.set_yticklabels(['%.3f' %flux for flux in list(flux_ticks)],fontsize=14.0)

	##Set some limits on the spextral axis - need to do it in logspace,
	##and set the edge gaps based on the data
	freq_delta = abs(np.log10(max(all_freqs))-np.log10(min(all_freqs)))/10.0
	flux_delta = abs(np.log10(max(all_fluxs))-np.log10(min(all_fluxs)))/10.0
	
	freq_max = 10**(np.log10(max(all_freqs))+freq_delta)
	freq_min = 10**(np.log10(min(all_freqs))-freq_delta)
	flux_max = 10**(np.log10(max(all_fluxs))+flux_delta)
	flux_min = 10**(np.log10(min(all_fluxs))-flux_delta)
	
	ax_spectral.set_xlim(freq_min,freq_max)
	ax_spectral.set_ylim(flux_min,flux_max)

	##Stick some grid stuff on the log log plot
	ten_steps = []
	for arr in [np.array([1e-5,2e-5,3e-5,4e-5,5e-5,6e-5,7e-5,8e-5,9e-5])*10**x for x in xrange(10)]:
		for i in list(arr): ten_steps.append(i)

	ax_spectral.set_xticks([step for step in ten_steps if step>freq_min and step<freq_max],minor=True)
	ax_spectral.set_yticks([step for step in ten_steps if step>flux_min and step<flux_max],minor=True)

	ax_spectral.xaxis.grid(True, which='minor',linestyle='dashed',alpha=0.1)
	ax_spectral.yaxis.grid(True, which='minor',linestyle='dashed',alpha=0.1)
	ax_spectral.set_xlabel(r'log$_{10}$(Frequency) (MHz)',fontsize=14)
	ax_spectral.set_ylabel(r'log$_{10}$(Flux) (Jy)',fontsize=14)

def create_plot(comp,accepted_inds,match_crit,dom_crit,outcome):
	
	###Split the information up as needed
	chunks = comp.split('START_COMP')
	all_info = chunks[0].split('\n')

	##FOR SOME REASON CAN'T DO BOTH OF THESE LINES IN THE SAME FOR LOOP?!?!?!
	for entry in all_info:   
		if entry=='': del all_info[all_info.index(entry)]
	for entry in all_info:
		if 'START' in entry: del all_info[all_info.index(entry)]

#	src_all = mkl.get_allinfo(all_info)

#	for i in xrange(len(src_all.cats)):
#		print "matched ['%s', %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f]" %(src_all.names[i],src_all.ras[i],src_all.rerrs[i],src_all.decs[i],src_all.derrs[i], src_all.freqs[i][0], src_all.fluxs[i][0],src_all.ferrs[i][0])

	matches = chunks[1].split('\n')
	del matches[0],matches[-2:]

	##See how many matches there are, and set up the number of plots needed
	num_matches = len(matches)
	if num_matches==1:
		width = 1
		height = 2
	else:
		width = int(num_matches**0.5)
		height = num_matches/width
		if num_matches%width==0:
			pass
		else:
			height+=1
	
	##Sets up a grid layout for the whole of the figure. We'll use half later on for
	##the individual plots
	gs = gridspec.GridSpec(height,2*width)
	
	##Need a big plot!
	fig = plt.figure(figsize = (18,11))

	##Find out the ra,dec of the base catalogue source
	info=all_info[0].split()
	ra_main = float(info[2])
	dec_main = float(info[4])
	
	main_dims = [0.16, 0.5, 0.29, 0.35]
	spec_dims = [0.16, 0.1, 0.29, 0.35]
	
	ax_main,ax_spectral,tr_fk5,wcs = make_left_plots(fig,main_dims,spec_dims,ra_main,dec_main)
	
	##Find the limits out to search area - have to do each edge individual,
	##because of the arcdistance projection malarky
	##Even at the same dec, 3 arcmins apart in RA doesn't translate to 3arcmin arcdist - projection
	##fun. Can use law of cosines on a sphere to work out appropriate delta RA. Use this to define plot
	##limits for a nice looking plot
	
	delta_RA = np.arccos((np.cos((2*closeness)*dr)-np.sin(dec_main*dr)**2)/np.cos(dec_main*dr)**2)/dr
	
	plot_lim = (2*closeness) + (0.1/60.0)
	ra_up_lim = ra_main + delta_RA + (0.1/60.0)
	ra_down_lim = ra_main - delta_RA - (0.1/60.0)
	dec_up_lim = dec_main + plot_lim
	dec_down_lim = dec_main - plot_lim
	
	##Plot the individual combination plots - do this first so the error bars go over
	##the top of the line plots
	spec_labels = []
	SIs = []
	for i in xrange(height):
		for j in range(width,2*width):
			try:
				ind = (i*width)+(j-width)
				match = matches[ind]
				ax = plt.subplot(gs[i,j])
				ax.set_xticklabels([])
				ax.set_yticklabels([])
				##TODO - if plot centred on or close to RA,Dec = 0,0 then going to get wrapping problems. Should be able to pull the need
				##for a wrap from ra_down_lim,ra_up_lim - one should be <0.0, or >360.0. Need to happen inside plot_ind
				prob,resid,spec_plot,params = plot_ind(match,ax,ind,ax_spectral,ra_down_lim,ra_up_lim,dec_down_lim,dec_up_lim,dom_crit,outcome)
				if spec_plot=='na':
					pass
				else:
					SIs.append([params[0],str(ind+1)])
					spec_labels.append(spec_plot)
			except IndexError:
				pass
	
	#===========================================================#
	##Plot the matching criteria information
	match1 = matches[0].split()
	src_g = mkl.get_srcg(match1)
	
	text_axes = fig.add_axes([0.45,0.5,0.125,0.35])
	text_axes.axis('off')

	##Plot the matching information
	props = dict(boxstyle='round', facecolor='w',lw='1.5')
	text_axes.text(0.5,0.5,'Match Criteria:\n%s\n\nDominace Test:\n%s\n\nOutcome:\n%s' %(match_crit,dom_crit,outcome),
		bbox=props,transform=text_axes.transAxes,verticalalignment='center',horizontalalignment='center',fontsize=16)
	
	all_fluxs = []
	
	##If no repeated catalogues to combine, skip
	if num_matches==0 or num_matches==1:
		pass
	##Otherwise, plot the combined fluxes
	else:
		##Calculate and plot the combined fluxes of the two sources, even if source one or two has been accepted
		##just as a guide
		src_all = mkl.get_allinfo(all_info)
		
		if accepted_inds=='Nope':
			pass
		else:
			comb_crit, ra_ws, rerr_ws, dec_ws, derr_ws, temp_freqs, comb_freqs, comb_fluxs, comb_ferrs, comb_fit, comb_jstat, comb_chi_red, combined_names, set_freqs, set_fluxs, set_fits, comb_sources = mkl.combine_flux(src_all,src_g,accepted_inds,'plot=yes',len(matches))
		
		##If the criteria sent the double to be combined, actually plot the fitted line
		if dom_crit == 'No dom. source':
			
			#for freq,flux in zip(set_freqs,set_fluxs):
				#ax_spectral.plot(freq,flux,linestyle='--',linewidth=1,color='r')
			
			split_colors = ['#AE70ED','#FFB60B','#62A9FF','#59DF00']
			for fit in set_fits:
				ind = set_fits.index(fit)
				#ax_spectral.plot(set_freqs[ind],set_fluxs[ind],linestyle='--',linewidth=1.0,color=split_colors[ind],alpha=0.7)
				split_p, = ax_spectral.plot(temp_freqs,np.exp(fit.params[1] + np.log(temp_freqs)*fit.params[0]),linestyle='-',linewidth=1.5,color=split_colors[ind])
				spec_labels.append(split_p)
				SIs.append([fit.params[0],'split %d' %(ind+1)])
			
			bright_colours = ['#FF6600','#33FF33','#FF47A3','#00ebb3']
			
			
			for freq in xrange(len(comb_freqs)):
				plot_errors('*',bright_colours[freq],comb_freqs[freq],comb_fluxs[freq],comb_ferrs[freq],'combo',9,ax_spectral)
			comb_p, = ax_spectral.plot(temp_freqs,np.exp(comb_fit.fittedvalues),linestyle='--',linewidth=1.5,color='k')
			spec_labels.append(comb_p)
			SIs.append([comb_fit.params[0],'comb'])
			
			##Send the combined fluxes to the all_fluxs so that ax_spectral is scaled appropriately
			for flux in comb_fluxs:
				all_fluxs.append(flux)
			
			#for pos in xrange(len(ra_ws)):
				#patch = plot_pos('*',bright_colours[pos],ra_ws[pos],dec_ws[pos],rerr_ws[pos],derr_ws[pos],combined_names[pos],14,ax_main,ax_main.get_transform("fk5"))

	##==============================================================

	##Fill the left hand plots with information goodness
	fill_left_plots(all_info,ra_main,dec_main,ax_main,ax_spectral,tr_fk5,wcs,all_fluxs,ra_down_lim,ra_up_lim,dec_down_lim,dec_up_lim,delta_RA)
	
		##If no repeated catalogues to combine, skip
	if num_matches==0 or num_matches==1:
		pass
	##Otherwise, plot the combined fluxes
	else:
		##Calculate and plot the combined fluxes of the two sources, even if source one or two has been accepted
		##just as a guide
		src_all = mkl.get_allinfo(all_info)
		bright_colours = ['#FF6600','#33FF33','#FF47A3','#00ebb3']
		if accepted_inds=='Nope':
			pass
		else:
			comb_crit, ra_ws, rerr_ws, dec_ws, derr_ws, temp_freqs, comb_freqs, comb_fluxs, comb_ferrs, comb_fit, comb_jstat, comb_chi_red, combined_names, set_freqs, set_fluxs, set_fits, comb_sources = mkl.combine_flux(src_all,src_g,accepted_inds,'plot=yes',len(matches))
		
			for pos in xrange(len(ra_ws)):
				patch = plot_pos('*',bright_colours[pos],ra_ws[pos],dec_ws[pos],0,0,combined_names[pos],14,ax_main,ax_main.get_transform("fk5")) #rerr_ws[pos],derr_ws[pos]
	
	fig.tight_layout()
	fig.subplots_adjust(bottom=0.1)
	fig.subplots_adjust(left=0.15)
	
	##Make room at the top of the plot for a legend for ax_main, make the legend
	fig.subplots_adjust(top=0.85)
	
	leg_labels = [r'$\alpha_{%s}$ = %.2f' %(SI[1],SI[0]) for SI in SIs]
	main_handles,main_labels = ax_main.get_legend_handles_labels()
	
	main_leg = fig.add_axes([0.05,0.87,0.9,0.05])
	main_leg.axis('off')
	main_leg.legend(main_handles,main_labels,loc='lower center',prop={'size':12},ncol=8) #,bbox_to_anchor=(0,1.02),
	
	spec_leg = fig.add_axes([0.45,0.1,0.125,0.35])
	spec_leg.axis('off')

	##Stop the legend from having so many entries that it goes off the plot	
	if len(spec_labels)>11:
		trim_labels = spec_labels[:10]
		trim_labels.append(spec_labels[-1])
		trim_legs = leg_labels[:10]
		trim_legs.append(leg_labels[-1])
		spec_leg.legend(trim_labels,trim_legs,loc='center',prop={'size':14},fancybox=True)
	else:
		spec_leg.legend(spec_labels,leg_labels,loc='center',prop={'size':14},fancybox=True)
	
	##Create an axes to contain patches for an ellipse legend
	patch_leg = fig.add_axes([0.015,0.1,0.06,0.75])
	patch_leg.set_xticks([])
	patch_leg.set_yticks([])
	patch_leg.set_xticklabels([])
	patch_leg.set_yticklabels([])
	
	##See what catalogues are present in the match
	present_cats = [cat for cat in set(src_g.cats) if cat!='-100000.0']
	##Scale accordingly
	increment = 1.0/(2+len(present_cats))
	ell_positions = np.arange(increment/2,1,increment)
	
	##Find the axes coord transform
	patch_trans = patch_leg.transAxes
	##Plot and name the resolution ellipse
	ell = patches.Ellipse((0.5,ell_positions[-2]),0.9,increment-0.05,angle=0,
		transform=patch_trans, linestyle='dashed',fc='none',lw=1.1,color='gray')
	patch_leg.add_patch(ell)
	patch_leg.text(0.5,ell_positions[-2],'Resolution\n+ error',
		transform=patch_trans,verticalalignment='center',horizontalalignment='center',fontsize=14)
	
	##Plot and name the search ellipse
	ell = patches.Ellipse((0.5,ell_positions[-1]),0.9,increment-0.05,angle=0,
		transform=patch_trans, linestyle='dashdot',fc='none',lw=1.1,color='k')
	patch_leg.add_patch(ell)
	patch_leg.text(0.5,ell_positions[-1],'Search\nradius',
		transform=patch_trans,verticalalignment='center',horizontalalignment='center',fontsize=14)
	
	##Use the same method as plot_all - for some reason was getting transform errors.
	##so do it separately here (sigh)
	for cat in present_cats:
		col_ind = matched_cats.index(cat)
		position_ind = present_cats.index(cat)
		patch_leg.errorbar(0.5,ell_positions[position_ind],0.01,0.075,marker=markers[col_ind],ms=8,mfc=marker_colours[col_ind],
			mec=marker_colours[col_ind],ecolor=marker_colours[col_ind],markeredgewidth=1,label='meh',linestyle='None',transform=patch_trans)
		
		ell = patches.Ellipse((0.5,ell_positions[position_ind]),0.9,increment-0.05,angle=0, transform=patch_trans,
			fc=ell_colours1[col_ind],color=ell_colours2[col_ind],alpha=alphas[col_ind])
		patch_leg.add_patch(ell)
		
		patch_leg.text(0.5,ell_positions[position_ind]-(increment/2-0.04),cat,
			transform=patch_trans,verticalalignment='center',horizontalalignment='center',fontsize=16)
	
	plt.show()
#	src_all = mkl.get_allinfo(all_info)
#	fig.savefig('bigchi_flatSI_%s_puma.png' %src_all.names[0],bbox_inches='tight',dpi=100)


def create_plot_fig(comp,accepted_inds,dom_crit,outcome):
	
	###Split the information up as needed
	chunks = comp.split('START_COMP')
	all_info = chunks[0].split('\n')
	
	##FOR SOME REASON CAN'T DO BOTH OF THESE LINES IN THE SAME FOR LOOP?!?!?!
	for entry in all_info:   
		if entry=='': del all_info[all_info.index(entry)]
	for entry in all_info:
		if 'START' in entry: del all_info[all_info.index(entry)]

	matches = chunks[1].split('\n')
	del matches[0],matches[-2:]

	##See how many matches there are, and set up the number of plots needed
	num_matches = len(matches)
	if num_matches==1:
		width = 1
		height = 2
	else:
		width = int(num_matches**0.5)
		height = num_matches/width
		if num_matches%width==0:
			pass
		else:
			height+=1
	
	##Sets up a grid layout for the whole of the figure. We'll use half later on for
	##the individual plots
	gs = gridspec.GridSpec(height,2*width)
	
	##Need a big plot!
	fig = plt.figure(figsize = (18,11))

	##Find out the ra,dec of the base catalogue source
	info=all_info[0].split()
	ra_main = float(info[2])
	dec_main = float(info[4])
	
	main_dims = [0.16, 0.5, 0.29, 0.35]
	spec_dims = [0.16, 0.1, 0.29, 0.35]
	
	ax_main,ax_spectral,tr_fk5,wcs = make_left_plots(fig,main_dims,spec_dims,ra_main,dec_main)
	
	##Find the limits out to search area - have to do each edge individual,
	##because of the arcdistance projection malarky
	##Even at the same dec, 3 arcmins apart in RA doesn't translate to 3arcmin arcdist - projection
	##fun. Can use law of cosines on a sphere to work out appropriate delta RA. Use this to define plot
	##limits for a nice looking plot
	
	delta_RA = np.arccos((np.cos((2*closeness)*dr)-np.sin(dec_main*dr)**2)/np.cos(dec_main*dr)**2)/dr
	
	plot_lim = (2*closeness) + (0.1/60.0)
	ra_up_lim = ra_main + delta_RA + (0.1/60.0)
	ra_down_lim = ra_main - delta_RA - (0.1/60.0)
	dec_up_lim = dec_main + plot_lim
	dec_down_lim = dec_main - plot_lim
	
	##Plot the individual combination plots - do this first so the error bars go over
	##the top of the line plots
	spec_labels = []
	SIs = []
	for i in xrange(height):
		for j in range(width,2*width):
			try:
				ind = (i*width)+(j-width)
				match = matches[ind]
				ax = plt.subplot(gs[i,j])
				ax.set_xticklabels([])
				ax.set_yticklabels([])
				##TODO - if plot centred on or close to RA,Dec = 0,0 then going to get wrapping problems. Should be able to pull the need
				##for a wrap from ra_down_lim,ra_up_lim - one should be <0.0, or >360.0. Need to happen inside plot_ind
				prob,resid,spec_plot,params = plot_ind(match,ax,ind,ax_spectral,ra_down_lim,ra_up_lim,dec_down_lim,dec_up_lim,dom_crit,outcome)
				if spec_plot=='na':
					pass
				else:
					SIs.append([params[0],str(ind+1)])
					spec_labels.append(spec_plot)
			except IndexError:
				pass
	
	#===========================================================#
	##Plot the matching criteria information
	match1 = matches[0].split()
	src_g = mkl.get_srcg(match1)
	
	text_axes = fig.add_axes([0.45,0.5,0.125,0.35])
	text_axes.axis('off')

	##Plot the matching information
	props = dict(boxstyle='round', facecolor='w',lw='1.5')
	text_axes.text(0.5,0.5,'Dominace Test:\n%s\n\nOutcome:\n%s' %(dom_crit,outcome),
		bbox=props,transform=text_axes.transAxes,verticalalignment='center',horizontalalignment='center',fontsize=16)
	
	all_fluxs = []
	
	##If no repeated catalogues to combine, skip
	if num_matches==0 or num_matches==1:
		pass
	##Otherwise, plot the combined fluxes
	else:
		##Calculate and plot the combined fluxes of the two sources, even if source one or two has been accepted
		##just as a guide
		src_all = mkl.get_allinfo(all_info)
		
		if accepted_inds=='Nope':
			pass
		else:
			comb_crit, ra_ws, rerr_ws, dec_ws, derr_ws, temp_freqs, comb_freqs, comb_fluxs, comb_ferrs, comb_fit, comb_jstat, comb_chi_red, combined_names, set_freqs, set_fluxs, set_fits, comb_sources = mkl.combine_flux(src_all,src_g,accepted_inds,'plot=yes',len(matches))
		
		##If the criteria sent the double to be combined, actually plot the fitted line
		if 'combine' in dom_crit or 'split' in dom_crit:
			
			#for freq,flux in zip(set_freqs,set_fluxs):
				#ax_spectral.plot(freq,flux,linestyle='--',linewidth=1,color='r')
			
			split_colors = ['#AE70ED','#FFB60B','#62A9FF','#59DF00']
			for fit in set_fits:
				ind = set_fits.index(fit)
				#ax_spectral.plot(set_freqs[ind],set_fluxs[ind],linestyle='--',linewidth=1.0,color=split_colors[ind],alpha=0.7)
				split_p, = ax_spectral.plot(temp_freqs,np.exp(fit.params[1] + np.log(temp_freqs)*fit.params[0]),linestyle='-',linewidth=1.5,color=split_colors[ind])
				spec_labels.append(split_p)
				SIs.append([fit.params[0],'split %d' %(ind+1)])
			
			bright_colours = ['#FF6600','#33FF33','#FF47A3','#00ebb3']
			
			for freq in xrange(len(comb_freqs)):
				plot_errors_comb('*',bright_colours[freq],comb_freqs[freq],comb_fluxs[freq],comb_ferrs[freq],'combo',16,ax_spectral)
			comb_p, = ax_spectral.plot(temp_freqs,np.exp(comb_fit.fittedvalues),linestyle='--',linewidth=1.5,color='k')
			spec_labels.append(comb_p)
			SIs.append([comb_fit.params[0],'comb'])
			
			##Send the combined fluxes to the all_fluxs so that ax_spectral is scaled appropriately
			for flux in comb_fluxs:
				all_fluxs.append(flux)
			
			#for pos in xrange(len(ra_ws)):
				#patch = plot_pos('*',bright_colours[pos],ra_ws[pos],dec_ws[pos],rerr_ws[pos],derr_ws[pos],combined_names[pos],14,ax_main,ax_main.get_transform("fk5"))

	##==============================================================

	##Fill the left hand plots with information goodness
	fill_left_plots(all_info,ra_main,dec_main,ax_main,ax_spectral,tr_fk5,wcs,all_fluxs,ra_down_lim,ra_up_lim,dec_down_lim,dec_up_lim,delta_RA)
	
		##If no repeated catalogues to combine, skip
	if num_matches==0 or num_matches==1:
		pass
	##Otherwise, plot the combined fluxes
	else:
		##Calculate and plot the combined fluxes of the two sources, even if source one or two has been accepted
		##just as a guide
		src_all = mkl.get_allinfo(all_info)
		bright_colours = ['#FF6600','#33FF33','#FF47A3','#00ebb3']
		if accepted_inds=='Nope':
			pass
		else:
			comb_crit, ra_ws, rerr_ws, dec_ws, derr_ws, temp_freqs, comb_freqs, comb_fluxs, comb_ferrs, comb_fit, comb_jstat, comb_chi_red, combined_names, set_freqs, set_fluxs, set_fits, comb_sources = mkl.combine_flux(src_all,src_g,accepted_inds,'plot=yes',len(matches))
		
			for pos in xrange(len(ra_ws)):
				patch = plot_pos('*',bright_colours[pos],ra_ws[pos],dec_ws[pos],0,0,combined_names[pos],14,ax_main,ax_main.get_transform("fk5")) #rerr_ws[pos],derr_ws[pos]
	
	fig.tight_layout()
	fig.subplots_adjust(bottom=0.1)
	fig.subplots_adjust(left=0.15)
	
	##Make room at the top of the plot for a legend for ax_main, make the legend
	fig.subplots_adjust(top=0.85)
	
	leg_labels = [r'$\alpha_{%s}$ = %.2f' %(SI[1],SI[0]) for SI in SIs]
	main_handles,main_labels = ax_main.get_legend_handles_labels()
	
	main_leg = fig.add_axes([0.05,0.87,0.9,0.05])
	main_leg.axis('off')
	main_leg.legend(main_handles,main_labels,loc='lower center',prop={'size':12},ncol=8) #,bbox_to_anchor=(0,1.02),
	
	spec_leg = fig.add_axes([0.45,0.1,0.125,0.35])
	spec_leg.axis('off')

	##Stop the legend from having so many entries that it goes off the plot	
	if len(spec_labels)>11:
		trim_labels = spec_labels[:10]
		trim_labels.append(spec_labels[-1])
		trim_legs = leg_labels[:10]
		trim_legs.append(leg_labels[-1])
		spec_leg.legend(trim_labels,trim_legs,loc='center',prop={'size':14},fancybox=True)
	else:
		spec_leg.legend(spec_labels,leg_labels,loc='center',prop={'size':14},fancybox=True)
	
	##Create an axes to contain patches for an ellipse legend
	patch_leg = fig.add_axes([0.015,0.1,0.06,0.75])
	patch_leg.set_xticks([])
	patch_leg.set_yticks([])
	patch_leg.set_xticklabels([])
	patch_leg.set_yticklabels([])
	
	##See what catalogues are present in the match
	present_cats = [cat for cat in set(src_g.cats) if cat!='-100000.0']
	##Scale accordingly
	increment = 1.0/(2+len(present_cats))
	ell_positions = np.arange(increment/2,1,increment)
	
	##Find the axes coord transform
	patch_trans = patch_leg.transAxes
	##Plot and name the resolution ellipse
	ell = patches.Ellipse((0.5,ell_positions[-2]),0.9,increment-0.05,angle=0,
		transform=patch_trans, linestyle='dashed',fc='none',lw=1.1,color='gray')
	patch_leg.add_patch(ell)
	patch_leg.text(0.5,ell_positions[-2],'Resolution\n+ error',
		transform=patch_trans,verticalalignment='center',horizontalalignment='center',fontsize=14)
	
	##Plot and name the search ellipse
	ell = patches.Ellipse((0.5,ell_positions[-1]),0.9,increment-0.05,angle=0,
		transform=patch_trans, linestyle='dashdot',fc='none',lw=1.1,color='k')
	patch_leg.add_patch(ell)
	patch_leg.text(0.5,ell_positions[-1],'Search\nradius',
		transform=patch_trans,verticalalignment='center',horizontalalignment='center',fontsize=14)
	
	##Use the same method as plot_all - for some reason was getting transform errors.
	##so do it separately here (sigh)
	for cat in present_cats:
		col_ind = matched_cats.index(cat)
		position_ind = present_cats.index(cat)
		patch_leg.errorbar(0.5,ell_positions[position_ind],0.01,0.075,marker=markers[col_ind],ms=8,mfc=marker_colours[col_ind],
			mec=marker_colours[col_ind],ecolor=marker_colours[col_ind],markeredgewidth=1,label='meh',linestyle='None',transform=patch_trans)
		
		ell = patches.Ellipse((0.5,ell_positions[position_ind]),0.9,increment-0.05,angle=0, transform=patch_trans,
			fc=ell_colours1[col_ind],color=ell_colours2[col_ind],alpha=alphas[col_ind])
		patch_leg.add_patch(ell)
		
		patch_leg.text(0.5,ell_positions[position_ind]-(increment/2-0.04),cat,
			transform=patch_trans,verticalalignment='center',horizontalalignment='center',fontsize=16)
	
	return fig


def create_plot_nocombo(comp,comp2,accepted_inds,num_matches,this_match,match_crit,dom_crit,outcome,tag,good_names):
	
	###Split the information up as needed
	chunks = comp.split('START_COMP')
	all_info = chunks[0].split('\n')

	##FOR SOME REASON CAN'T DO BOTH OF THESE LINES IN THE SAME FOR LOOP?!?!?!
	for entry in all_info:   
		if entry=='': del all_info[all_info.index(entry)]
	for entry in all_info:
		if 'START' in entry: del all_info[all_info.index(entry)]
		
	chunks2 = comp2.split('START_COMP')
	all_info2 = chunks2[0].split('\n')

	##FOR SOME REASON CAN'T DO BOTH OF THESE LINES IN THE SAME FOR LOOP?!?!?!
	for entry in all_info2:   
		if entry=='': del all_info2[all_info2.index(entry)]
	for entry in all_info2:
		if 'START' in entry: del all_info2[all_info2.index(entry)]

#	src_all = mkl.get_allinfo(all_info)

#	for i in xrange(len(src_all.cats)):
#		print "matched ['%s', %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f]" %(src_all.names[i],src_all.ras[i],src_all.rerrs[i],src_all.decs[i],src_all.derrs[i], src_all.freqs[i][0], src_all.fluxs[i][0],src_all.ferrs[i][0])

	matches = chunks[1].split('\n')
	del matches[0],matches[-2:]

	###See how many matches there are, and set up the number of plots needed
	#num_matches = len(matches)
	#if num_matches==1:
		#width = 1
		#height = 2
	#else:
		#width = int(num_matches**0.5)
		#height = num_matches/width
		#if num_matches%width==0:
			#pass
		#else:
			#height+=1
	
	###Sets up a grid layout for the whole of the figure. We'll use half later on for
	###the individual plots
	#gs = gridspec.GridSpec(height,2*width)
	
	##Need a big plot!
	fig = plt.figure(figsize = (14,11))

	##Find out the ra,dec of the base catalogue source
	info=all_info[0].split()
	ra_main = float(info[2])
	dec_main = float(info[4])
	
	main_dims = [0.18, 0.5, 0.32, 0.35]
	spec_dims = [0.18, 0.1, 0.32, 0.35]
	
	ax_main,ax_spectral,tr_fk5,wcs = make_left_plots(fig,main_dims,spec_dims,ra_main,dec_main)
	
	
	main_dims_right = [0.55, 0.5, 0.32, 0.35]
	spec_dims_right = [0.55, 0.1, 0.32, 0.35]
	
	ax_main_right,ax_spectral_right,tr_fk5_right,wcs_right = make_right_plots(fig,main_dims_right,spec_dims_right,ra_main,dec_main)
	
	##Find the limits out to search area - have to do each edge individual,
	##because of the arcdistance projection malarky
	##Even at the same dec, 3 arcmins apart in RA doesn't translate to 3arcmin arcdist - projection
	##fun. Can use law of cosines on a sphere to work out appropriate delta RA. Use this to define plot
	##limits for a nice looking plot
	
	delta_RA = np.arccos((np.cos((2*closeness)*dr)-np.sin(dec_main*dr)**2)/np.cos(dec_main*dr)**2)/dr
	
	plot_lim = (2*closeness) + (0.1/60.0)
	ra_up_lim = ra_main + delta_RA + (0.1/60.0)
	ra_down_lim = ra_main - delta_RA - (0.1/60.0)
	dec_up_lim = dec_main + plot_lim
	dec_down_lim = dec_main - plot_lim
	
	###Plot the individual combination plots - do this first so the error bars go over
	###the top of the line plots
	spec_labels = []
	SIs = []
	#for i in xrange(height):
		#for j in range(width,2*width):
			#try:
				#ind = (i*width)+(j-width)
				#match = matches[ind]
				#ax = plt.subplot(gs[i,j])
				#ax.set_xticklabels([])
				#ax.set_yticklabels([])
				###TODO - if plot centred on or close to RA,Dec = 0,0 then going to get wrapping problems. Should be able to pull the need
				###for a wrap from ra_down_lim,ra_up_lim - one should be <0.0, or >360.0. Need to happen inside plot_ind
				#prob,resid,spec_plot,params = plot_ind(match,ax,ind,ax_spectral,ra_down_lim,ra_up_lim,dec_down_lim,dec_up_lim,dom_crit,outcome)
				#if spec_plot=='na':
					#pass
				#else:
					#SIs.append([params[0],str(ind+1)])
					#spec_labels.append(spec_plot)
			#except IndexError:
				#pass
	
	#===========================================================#
	##Plot the matching criteria information
	match1 = matches[0].split()
	src_g = mkl.get_srcg(match1)
	
	#text_axes = fig.add_axes([0.45,0.5,0.125,0.35])
	#text_axes.axis('off')

	##Plot the matching information
	#props = dict(boxstyle='round', facecolor='w',lw='1.5')
	#text_axes.text(0.5,0.5,'Original\nMatch Criteria:\n%s\n\nDominace Test:\n%s\n\nOutcome:\n%s' %(match_crit,dom_crit,outcome),
		#bbox=props,transform=text_axes.transAxes,verticalalignment='center',horizontalalignment='center',fontsize=16)
	
	all_fluxs = []
	
	##If no repeated catalogues to combine, skip
	if num_matches==0:
		pass
	elif num_matches==1:
		for match in matches:
			if 'END_COMP' in match:
				pass
			else:
				ind = matches.index(match)
				prob,resid,spec_plot,params = plot_ind(match,ind,ax_spectral,dom_crit,outcome)
				#plot_ind(match,ax,ind,ax_spectral,ra_down_lim,ra_up_lim,dec_down_lim,dec_up_lim,dom_crit,outcome)
				if spec_plot=='na':
					pass
				else:
					SIs.append([params[0],str(ind+1)])
					spec_labels.append(spec_plot)
		
	##Otherwise, plot the combined fluxes
	else:
		##Calculate and plot the combined fluxes of the two sources, even if source one or two has been accepted
		##just as a guide
		src_all = mkl.get_allinfo(all_info)
		
		if accepted_inds=='Nope':
			pass
		else:
			comb_crit, ra_ws, rerr_ws, dec_ws, derr_ws, temp_freqs, comb_freqs, comb_fluxs, comb_ferrs, comb_fit, comb_jstat, comb_chi_red, combined_names, set_freqs, set_fluxs, set_fits, comb_sources = mkl.combine_flux(src_all,src_g,accepted_inds,'plot=yes',len(matches))
		
		##If the criteria sent the double to be combined, actually plot the fitted line
		#if dom_crit == 'No dom. source':
			
			#for freq,flux in zip(set_freqs,set_fluxs):
				#ax_spectral.plot(freq,flux,linestyle='--',linewidth=1,color='r')
			
			#split_colors = ['#AE70ED','#FFB60B','#62A9FF','#59DF00']
			#for fit in set_fits:
				#ind = set_fits.index(fit)
				##ax_spectral.plot(set_freqs[ind],set_fluxs[ind],linestyle='--',linewidth=1.0,color=split_colors[ind],alpha=0.7)
				#split_p, = ax_spectral.plot(temp_freqs,np.exp(fit.params[1] + np.log(temp_freqs)*fit.params[0]),linestyle='-',linewidth=1.5,color=split_colors[ind])
				#spec_labels.append(split_p)
				#SIs.append([fit.params[0],'split %d' %(ind+1)])
			
			bright_colours = ['#FF6600','#33FF33','#FF47A3','#00ebb3']
			
			for freq in xrange(len(comb_freqs)):
				plot_errors_comb('*',bright_colours[freq],comb_freqs[freq],comb_fluxs[freq],comb_ferrs[freq],'combo',16,ax_spectral)
			comb_p, = ax_spectral.plot(temp_freqs,np.exp(comb_fit.fittedvalues),linestyle='--',linewidth=1.5,color='k')
			spec_labels.append(comb_p)
			SIs.append([comb_fit.params[0],'comb'])
			
			##Send the combined fluxes to the all_fluxs so that ax_spectral is scaled appropriately
			for flux in comb_fluxs:
				all_fluxs.append(flux)
			
			#for pos in xrange(len(ra_ws)):
				#patch = plot_pos('*',bright_colours[pos],ra_ws[pos],dec_ws[pos],rerr_ws[pos],derr_ws[pos],combined_names[pos],14,ax_main,ax_main.get_transform("fk5"))

	##==============================================================
	matches2 = chunks[1].split('\n')
	del matches2[0],matches2[-2:]
	match12 = matches2[0].split()
	
	src_all2 = get_allinfo(all_info2,good_names)
	src_g2 = get_srcg(src_all2,match12,good_names)
	
	if this_match != 'Nope':
		pass
	else:
		comb_crit, ra_ws, rerr_ws, dec_ws, derr_ws, temp_freqs, comb_freqs, comb_fluxs, comb_ferrs, comb_fit, comb_jstat, comb_chi_red, combined_names, set_freqs, set_fluxs, set_fits, comb_sources = mkl.combine_flux(src_all2,src_g2,accepted_inds,'plot=yes',len(matches))
	
		bright_colours = ['#FF6600','#33FF33','#FF47A3','#00ebb3']
		
		for freq in xrange(len(comb_freqs)):
			plot_errors_comb('*',bright_colours[freq],comb_freqs[freq],comb_fluxs[freq],comb_ferrs[freq],'combo',16,ax_spectral_right)
		comb_p, = ax_spectral_right.plot(temp_freqs,np.exp(comb_fit.fittedvalues),linestyle='--',linewidth=1.5,color='m')
		spec_labels.append(comb_p)
		SIs.append([comb_fit.params[0],'comb_{fix}'])
		
		##Send the combined fluxes to the all_fluxs so that ax_spectral is scaled appropriately
		for flux in comb_fluxs:
			all_fluxs.append(flux)
	

	##Fill the left hand plots with information goodness
	fill_left_plots(all_info,ra_main,dec_main,ax_main,ax_spectral,tr_fk5,wcs,all_fluxs,ra_down_lim,ra_up_lim,dec_down_lim,dec_up_lim,delta_RA)
	fill_right_plots(all_info2,ra_main,dec_main,ax_main_right,ax_spectral_right,tr_fk5_right,wcs_right,all_fluxs,ra_down_lim,ra_up_lim,dec_down_lim,dec_up_lim,delta_RA,good_names)
	
		##If no repeated catalogues to combine, skip
	if num_matches==0 or num_matches==1:
		pass
	##Otherwise, plot the combined fluxes
	else:
		##Calculate and plot the combined fluxes of the two sources, even if source one or two has been accepted
		##just as a guide
		src_all = mkl.get_allinfo(all_info)
		bright_colours = ['#FF6600','#33FF33','#FF47A3','#00ebb3']
		if accepted_inds=='Nope':
			pass
		else:
			comb_crit, ra_ws, rerr_ws, dec_ws, derr_ws, temp_freqs, comb_freqs, comb_fluxs, comb_ferrs, comb_fit, comb_jstat, comb_chi_red, combined_names, set_freqs, set_fluxs, set_fits, comb_sources = mkl.combine_flux(src_all,src_g,accepted_inds,'plot=yes',len(matches))
		
			for pos in xrange(len(ra_ws)):
				patch = plot_pos_comb('*',bright_colours[pos],ra_ws[pos],dec_ws[pos],0,0,combined_names[pos],16,ax_main,ax_main.get_transform("fk5")) #rerr_ws[pos],derr_ws[pos]
				
	if this_match != 'Nope':
		accepted_match = matches2[int(this_match)-1].split()
		
		this_srcg = mkl.get_srcg(accepted_match)
		
		#fluxs = [flux for flux in this_srcg.fluxs if flux != -100000.0]
		#freqs = [freq for freq in this_srcg.freqs if freq != -100000.0]
		#ferrs = [ferr for ferr in this_srcg.ferrs if ferr != -100000.0]
		
		temp_freqs = []
		temp_fluxs = []
		temp_ferrs = []
	
		##Get the sources out of the array in list format (which is used later when making the sources
		##to add to the final table)
		
		for i in xrange(len(this_srcg.freqs)):
			for j in xrange(len(this_srcg.freqs[i])):
				if this_srcg.fluxs[i][j] == -100000.0 or np.isnan(this_srcg.fluxs[i][j])==True:
					pass
				else:
					temp_freqs.append(this_srcg.freqs[i][j])
		for i in xrange(len(this_srcg.freqs)):
			for j in xrange(len(this_srcg.freqs[i])):
				if this_srcg.fluxs[i][j] == -100000.0 or np.isnan(this_srcg.fluxs[i][j])==True:
					pass
				else:
					temp_fluxs.append(this_srcg.fluxs[i][j])
		for i in xrange(len(this_srcg.freqs)):
			for j in xrange(len(this_srcg.freqs[i])):
				if this_srcg.fluxs[i][j] == -100000.0 or np.isnan(this_srcg.fluxs[i][j])==True:
					pass
				else:
					temp_ferrs.append(this_srcg.ferrs[i][j]/this_srcg.fluxs[i][j])
		
		fluxs, freqs, ferrs = temp_fluxs, temp_freqs, temp_ferrs
		
		log_fluxs = np.log([flux for (freq,flux) in sorted(zip(freqs,fluxs),key=lambda pair: pair[0])])
		sorted_ferrs = np.array([ferr for (freq,ferr) in sorted(zip(freqs,ferrs),key=lambda pair: pair[0])])
		log_freqs = np.log(sorted(freqs))
		ferrs = np.array(ferrs)
		
		lin_fit,jstat,bse,chi_red = mkl.fit_line(log_freqs,log_fluxs,sorted_ferrs)
		#jstat_resids,params,bses,chi_resids = mkl.calculate_resids([accepted_match])
		
		spec_plot, = ax_spectral_right.plot(np.exp(log_freqs),np.exp(lin_fit.fittedvalues),linestyle='-',linewidth=1,alpha=0.7)
		spec_labels.append(spec_plot)
		SIs.append([lin_fit.params[0],'{fix}'])
		
	else:
		bright_colours = ['#FF6600','#33FF33','#FF47A3','#00ebb3']
		if accepted_inds=='Nope':
			pass
		else:
			comb_crit, ra_ws, rerr_ws, dec_ws, derr_ws, temp_freqs, comb_freqs, comb_fluxs, comb_ferrs, comb_fit, comb_jstat, comb_chi_red, combined_names, set_freqs, set_fluxs, set_fits, comb_sources = mkl.combine_flux(src_all2,src_g2,accepted_inds,'plot=yes',len(matches))
		
			for pos in xrange(len(ra_ws)):
				patch = plot_pos_comb('*',bright_colours[pos],ra_ws[pos],dec_ws[pos],0,0,combined_names[pos],16,ax_main_right,ax_main_right.get_transform("fk5")) #rerr_ws[pos],derr_ws[pos]
	
	#fig.tight_layout()
	fig.subplots_adjust(bottom=0.1)
	fig.subplots_adjust(left=0.15)
	
	##Make room at the top of the plot for a legend for ax_main, make the legend
	fig.subplots_adjust(top=0.85)
	
	leg_labels = [r'$\alpha_{%s}$ = %.2f' %(SI[1],SI[0]) for SI in SIs]
	main_handles,main_labels = ax_main.get_legend_handles_labels()

	main_right_handles,main_right_labels = ax_main_right.get_legend_handles_labels()
	
	for label in main_right_labels:
		if label not in main_labels and 'ombine' not in label:
			main_labels.append(label)
			main_handles.append(main_right_handles[main_right_labels.index(label)])
	
	
	main_leg = fig.add_axes([0.05,0.87,0.9,0.05])
	main_leg.axis('off')
	main_leg.legend(main_handles,main_labels,loc='lower center',prop={'size':12},ncol=8) #,bbox_to_anchor=(0,1.02),
	
	spec_leg = fig.add_axes([0.48,0.15,0.125,0.35])
	spec_leg.axis('off')

	##Stop the legend from having so many entries that it goes off the plot	
	if len(spec_labels)>11:
		trim_labels = spec_labels[:10]
		trim_labels.append(spec_labels[-1])
		trim_legs = leg_labels[:10]
		trim_legs.append(leg_labels[-1])
		spec_leg.legend(trim_labels,trim_legs,loc='center',prop={'size':14},fancybox=True)
	else:
		spec_leg.legend(spec_labels,leg_labels,loc='center',prop={'size':14},fancybox=True)
	
	##Create an axes to contain patches for an ellipse legend
	patch_leg = fig.add_axes([0.015,0.1,0.06,0.75])
	patch_leg.set_xticks([])
	patch_leg.set_yticks([])
	patch_leg.set_xticklabels([])
	patch_leg.set_yticklabels([])
	
	##See what catalogues are present in the match
	present_cats = [cat for cat in set(src_g.cats) if cat!='-100000.0']
	##Scale accordingly
	increment = 1.0/(2+len(present_cats))
	ell_positions = np.arange(increment/2,1,increment)
	
	##Find the axes coord transform
	patch_trans = patch_leg.transAxes
	##Plot and name the resolution ellipse
	ell = patches.Ellipse((0.5,ell_positions[-2]),0.9,increment-0.05,angle=0,
		transform=patch_trans, linestyle='dashed',fc='none',lw=1.1,color='gray')
	patch_leg.add_patch(ell)
	patch_leg.text(0.5,ell_positions[-2],'Resolution\n+ error',
		transform=patch_trans,verticalalignment='center',horizontalalignment='center',fontsize=14)
	
	##Plot and name the search ellipse
	ell = patches.Ellipse((0.5,ell_positions[-1]),0.9,increment-0.05,angle=0,
		transform=patch_trans, linestyle='dashdot',fc='none',lw=1.1,color='k')
	patch_leg.add_patch(ell)
	patch_leg.text(0.5,ell_positions[-1],'Search\nradius',
		transform=patch_trans,verticalalignment='center',horizontalalignment='center',fontsize=14)
	
	##Use the same method as plot_all - for some reason was getting transform errors.
	##so do it separately here (sigh)
	for cat in present_cats:
		col_ind = matched_cats.index(cat)
		position_ind = present_cats.index(cat)
		patch_leg.errorbar(0.5,ell_positions[position_ind],0.01,0.075,marker=markers[col_ind],ms=8,mfc=marker_colours[col_ind],
			mec=marker_colours[col_ind],ecolor=marker_colours[col_ind],markeredgewidth=1,label='meh',linestyle='None',transform=patch_trans)
		
		ell = patches.Ellipse((0.5,ell_positions[position_ind]),0.9,increment-0.05,angle=0, transform=patch_trans,
			fc=ell_colours1[col_ind],color=ell_colours2[col_ind],alpha=alphas[col_ind])
		patch_leg.add_patch(ell)
		
		if cat == 'nvss': cat = 'NVSS'
		elif cat == 'vlssr': cat = 'VLSSr'
		elif cat == 'comp_v11': cat = 'KGS'
		
		patch_leg.text(0.5,ell_positions[position_ind]-(increment/2-0.04),cat,
			transform=patch_trans,verticalalignment='center',horizontalalignment='center',fontsize=16)
	
	#plt.show()
	src_all = mkl.get_allinfo(all_info)
	
	
	
	fig.savefig('%s_edit.png' %(src_all.names[0]),bbox_inches='tight',dpi=100)
	#fig.savefig('paper_edit_top.pdf',bbox_inches='tight',dpi=200)
	
	#subprocess.call('pdftops -eps paper_edit_top.pdf paper_edit_top.eps',shell=True)
	#subprocess.call('rm paper_edit_top.pdf',shell=True)
	
	
	plt.close()
	
def combine_flux_nobase(src_all,src_g,accepted_inds,plot,num_matches):
	'''Takes a src_group() class that contains all group info. Indentifies which catalogue is repeated, combines
	the fluxes and fits a new line. Returns the reduced frequency list, the combined flux and flux error 
	arrays, as well as the line_fit object and residuals'''
	
	##Find repeated catalogues
	repeated_cats = set([src_all.cats[ind] for ind in accepted_inds if src_all.cats.count(src_all.cats[ind])>1])
	##This won't neccesarily be in the that the cats appear in src_all.cats so reorder
	repeat_indexs = [src_all.cats.index(cat) for cat in repeated_cats]
	repeated_cats = [cat for ind,cat in sorted(zip(repeat_indexs,repeated_cats),key=lambda pair: pair[0])]
	
	##These are used to test the combined spectrum
	#temp_freqs = [src_all.freqs[i] for i in xrange(len(src_all.freqs)) if src_all.cats[i] not in repeated_cats]
	#temp_fluxs = [src_all.fluxs[i] for i in xrange(len(src_all.fluxs)) if src_all.cats[i] not in repeated_cats]
	#temp_ferrs = [src_all.ferrs[i] for i in xrange(len(src_all.ferrs)) if src_all.cats[i] not in repeated_cats]
	
	temp_freqs = [src_all.freqs[i] for i in xrange(len(src_all.freqs)) if src_all.cats[i] not in repeated_cats and src_all.cats.index(src_all.cats[i])!=0]
	temp_fluxs = [src_all.fluxs[i] for i in xrange(len(src_all.fluxs)) if src_all.cats[i] not in repeated_cats and src_all.cats.index(src_all.cats[i])!=0]
	temp_ferrs = [src_all.ferrs[i] for i in xrange(len(src_all.ferrs)) if src_all.cats[i] not in repeated_cats and src_all.cats.index(src_all.cats[i])!=0]
	
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
		
		##TEST TO SEE IF THE REPEATED SOURCES ARE RESOLVED (BY A GIVEN RESOLUTION THRESHOLD)
		##Need to do the split test here even if not propagating to the final catalogue
		##-------------------------------------------------------------------------------------------------------
		
		##Give it a value even if not splitting
		if split==0:
			dist_test = 0.020833333 ##1.25 arcmin
		else:
			dist_test = dec_to_deg(split)
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
		temp_fluxs.append(comb_flux)
		temp_ferrs.append(comb_ferr)
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
						if arcdist(src_all.ras[src],src_all.ras[other_src],src_all.decs[src],src_all.decs[other_src])<dist_test:
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
			src_g.fluxs[srcg_ind] = comb_fluxs[i]
			src_g.ferrs[srcg_ind] = comb_ferrs[i]
		#srg_g.freqs = temp_freqs
		src_g.SI = float(comb_fit.params[0])
		src_g.intercept = comb_fit.params[1]
		src_g.SI_err = comb_bse[0]
		src_g.intercept_err = comb_bse[1]
	
		##If good fit, report that in the final stats object
		if comb_chi_red<=2:
			src_g.low_resids = 0
		else:
			src_g.low_resids = 1
			
		dom_crit = 'Accepted -\ncombined'
		
		if split != 0:
			##If a split source, create however many new sources are needed
			if big_sep=='yes':
			##Test to see if any of the new components fail a spec test
			##If so, flag out for eyeballing
				resid_tests = []
				for jstat,red in zip(set_jstat,set_red):
					if jstat<=jstat_thresh or red<=chi_thresh:
						resid_tests.append('yes')
					else:
						resid_tests.append('no')
				if 'no' in resid_tests:
					split_flag = 'split breaks it'
					dom_crit = 'Rejected -\nsplit'
				else:
					dom_crit = 'Accepted -\nsplit'
					split_sources = []
					for set_ind in xrange(len(set_cats)):
						new_g = copy.deepcopy(src_g)
						##We need to put the sources in the same order as the src_g, so it gets
						##put in to the final table in the right order. The position info for
						##everything is in the src_all info, and the order of the sources is in 
						##the src_g. There are also blank entries in src_g, so just make a copy
						##and insert the correct sources in to it.
						for src in xrange(len(set_cats[set_ind])):
							order_ind = src_g.cats.index(set_cats[set_ind][src])
							info_ind = src_all.names.index(set_names[set_ind][src])
							new_g.freqs[order_ind] = np.array([set_freqs[set_ind][src]])
							new_g.fluxs[order_ind] = np.array([set_fluxs[set_ind][src]])
							new_g.ferrs[order_ind] = np.array([set_ferrs[set_ind][src]])
							new_g.names[order_ind] = src_all.names[info_ind]
							new_g.ras[order_ind] = src_all.ras[info_ind]
							new_g.rerrs[order_ind] = src_all.rerrs[info_ind]
							new_g.decs[order_ind] = src_all.decs[info_ind]
							new_g.derrs[order_ind] = src_all.derrs[info_ind]
							new_g.minors[order_ind] = src_all.minors[info_ind]
							new_g.majors[order_ind] = src_all.majors[info_ind]
							new_g.PAs[order_ind] = src_all.PAs[info_ind]
							new_g.SI = set_fits[set_ind].params[0]
							new_g.SI_err = set_bse[set_ind][0]
							new_g.intercept = set_fits[set_ind].params[1]
							new_g.intercept_err = set_bse[set_ind][1]
							
						split_sources.append(new_g)
			
		##If plotting need a bunch specific info
		if plot=='plot=yes':
			##Get the freqs, fluxes and ferrs out of array for and in to a simple list to be plotted
			comb_freqs_sing = []
			comb_fluxs_sing = []
			comb_ferrs_sing = []
			for freqs in comb_freqs:
				for freq in freqs: comb_freqs_sing.append(freq)
			for fluxs in comb_fluxs:
				for flux in fluxs: comb_fluxs_sing.append(flux)
			for ferrs in comb_ferrs:
				for ferr in ferrs: comb_ferrs_sing.append(ferr)
				
			if dom_crit == 'Rejected -\nsplit':
				return dom_crit, ra_ws, rerr_ws, dec_ws, derr_ws, np.exp(log_temp_freqs), comb_freqs_sing, comb_fluxs_sing, comb_ferrs_sing, comb_fit, comb_jstat, comb_chi_red, combined_names, set_freqs, set_fluxs, set_fits,[src_g] 
			else:
				if 'combined' in dom_crit:
					return dom_crit, ra_ws, rerr_ws, dec_ws, derr_ws, np.exp(log_temp_freqs), comb_freqs_sing, comb_fluxs_sing, comb_ferrs_sing, comb_fit, comb_jstat, comb_chi_red, combined_names, set_freqs, set_fluxs, set_fits,[src_g] 
				else:
					return dom_crit, ra_ws, rerr_ws, dec_ws, derr_ws, np.exp(log_temp_freqs), comb_freqs_sing, comb_fluxs_sing, comb_ferrs_sing, comb_fit, comb_jstat, comb_chi_red, combined_names, set_freqs, set_fluxs, set_fits,split_sources
		else:
			if dom_crit == 'Rejected -\nsplit':
				##It's been failed by split, but had passed by combine
				return dom_crit, 'nyope', comb_jstat, comb_chi_red
			else:
				##It's passed the combine
				if 'combined' in dom_crit:
					##Return just the combined source
					return dom_crit, [src_g], comb_jstat, comb_chi_red
				else:
					##Return all of the components
					return dom_crit, split_sources, comb_jstat, comb_chi_red
	
	##If it fails, still return all the info to the plot so we can see what's going on
	else:
		if split==False:
			big_sep='no'
		
		if big_sep=='yes':
			dom_crit = 'To eyeball\n(split?)'
		else:
			dom_crit = 'To eyeball'
		
		if plot=='plot=yes':
			comb_freqs_sing = []
			comb_fluxs_sing = []
			comb_ferrs_sing = []
			for freqs in comb_freqs:
				for freq in freqs: comb_freqs_sing.append(freq)
			for fluxs in comb_fluxs:
				for flux in fluxs: comb_fluxs_sing.append(flux)
			for ferrs in comb_ferrs:
				for ferr in ferrs: comb_ferrs_sing.append(ferr)
				
			########################################################################################
				
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
				src_g.fluxs[srcg_ind] = comb_fluxs[i]
				src_g.ferrs[srcg_ind] = comb_ferrs[i]
			#srg_g.freqs = temp_freqs
			src_g.SI = float(comb_fit.params[0])
			src_g.intercept = comb_fit.params[1]
			src_g.SI_err = comb_bse[0]
			src_g.intercept_err = comb_bse[1]
		
			##If good fit, report that in the final stats object
			if comb_chi_red<=2:
				src_g.low_resids = 0
			else:
				src_g.low_resids = 1
						
			########################################################################################
				
			return dom_crit, ra_ws, rerr_ws, dec_ws, derr_ws, np.exp(log_temp_freqs), comb_freqs_sing, comb_fluxs_sing, comb_ferrs_sing, comb_fit, comb_jstat, comb_chi_red, ["(combined-%s)" %cat for cat in repeated_cats], set_freqs ,set_fluxs, set_fits,[src_g]
		else:
			return dom_crit, 'nyope', comb_jstat, comb_chi_red

def create_plot_nocombo_spilt(comp,comp2,accepted_inds,num_matches,match_crit,dom_crit,outcome,tag):
	
	###Split the information up as needed
	chunks = comp.split('START_COMP')
	all_info = chunks[0].split('\n')

	##FOR SOME REASON CAN'T DO BOTH OF THESE LINES IN THE SAME FOR LOOP?!?!?!
	for entry in all_info:   
		if entry=='': del all_info[all_info.index(entry)]
	for entry in all_info:
		if 'START' in entry: del all_info[all_info.index(entry)]
		
	chunks2 = comp2.split('START_COMP')
	all_info2 = chunks2[0].split('\n')

	##FOR SOME REASON CAN'T DO BOTH OF THESE LINES IN THE SAME FOR LOOP?!?!?!
	for entry in all_info2:   
		if entry=='': del all_info2[all_info2.index(entry)]
	for entry in all_info2:
		if 'START' in entry: del all_info2[all_info2.index(entry)]

#	src_all = mkl.get_allinfo(all_info)

#	for i in xrange(len(src_all.cats)):
#		print "matched ['%s', %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f]" %(src_all.names[i],src_all.ras[i],src_all.rerrs[i],src_all.decs[i],src_all.derrs[i], src_all.freqs[i][0], src_all.fluxs[i][0],src_all.ferrs[i][0])

	matches = chunks[1].split('\n')
	del matches[0],matches[-2:]

	###See how many matches there are, and set up the number of plots needed
	#num_matches = len(matches)
	#if num_matches==1:
		#width = 1
		#height = 2
	#else:
		#width = int(num_matches**0.5)
		#height = num_matches/width
		#if num_matches%width==0:
			#pass
		#else:
			#height+=1
	
	###Sets up a grid layout for the whole of the figure. We'll use half later on for
	###the individual plots
	#gs = gridspec.GridSpec(height,2*width)
	
	##Need a big plot!
	fig = plt.figure(figsize = (14,11))

	##Find out the ra,dec of the base catalogue source
	info=all_info[0].split()
	ra_main = float(info[2])
	dec_main = float(info[4])
	
	main_dims = [0.18, 0.5, 0.32, 0.35]
	spec_dims = [0.18, 0.1, 0.32, 0.35]
	
	ax_main,ax_spectral,tr_fk5,wcs = make_left_plots(fig,main_dims,spec_dims,ra_main,dec_main)
	
	
	main_dims_right = [0.55, 0.5, 0.32, 0.35]
	spec_dims_right = [0.55, 0.1, 0.32, 0.35]
	
	ax_main_right,ax_spectral_right,tr_fk5_right,wcs_right = make_right_plots(fig,main_dims_right,spec_dims_right,ra_main,dec_main)
	
	##Find the limits out to search area - have to do each edge individual,
	##because of the arcdistance projection malarky
	##Even at the same dec, 3 arcmins apart in RA doesn't translate to 3arcmin arcdist - projection
	##fun. Can use law of cosines on a sphere to work out appropriate delta RA. Use this to define plot
	##limits for a nice looking plot
	
	delta_RA = np.arccos((np.cos((3*closeness)*dr)-np.sin(dec_main*dr)**2)/np.cos(dec_main*dr)**2)/dr
	
	plot_lim = (3*closeness) + (0.1/60.0)
	ra_up_lim = ra_main + delta_RA + (0.1/60.0)
	ra_down_lim = ra_main - delta_RA - (0.1/60.0)
	dec_up_lim = dec_main + plot_lim
	dec_down_lim = dec_main - plot_lim
	
	###Plot the individual combination plots - do this first so the error bars go over
	###the top of the line plots
	spec_labels = []
	SIs = []
	#for i in xrange(height):
		#for j in range(width,2*width):
			#try:
				#ind = (i*width)+(j-width)
				#match = matches[ind]
				#ax = plt.subplot(gs[i,j])
				#ax.set_xticklabels([])
				#ax.set_yticklabels([])
				###TODO - if plot centred on or close to RA,Dec = 0,0 then going to get wrapping problems. Should be able to pull the need
				###for a wrap from ra_down_lim,ra_up_lim - one should be <0.0, or >360.0. Need to happen inside plot_ind
				#prob,resid,spec_plot,params = plot_ind(match,ax,ind,ax_spectral,ra_down_lim,ra_up_lim,dec_down_lim,dec_up_lim,dom_crit,outcome)
				#if spec_plot=='na':
					#pass
				#else:
					#SIs.append([params[0],str(ind+1)])
					#spec_labels.append(spec_plot)
			#except IndexError:
				#pass
	
	#===========================================================#
	##Plot the matching criteria information
	match1 = matches[0].split()
	src_g = mkl.get_srcg(match1)
	
	#text_axes = fig.add_axes([0.45,0.5,0.125,0.35])
	#text_axes.axis('off')

	##Plot the matching information
	#props = dict(boxstyle='round', facecolor='w',lw='1.5')
	#text_axes.text(0.5,0.5,'Original\nMatch Criteria:\n%s\n\nDominace Test:\n%s\n\nOutcome:\n%s' %(match_crit,dom_crit,outcome),
		#bbox=props,transform=text_axes.transAxes,verticalalignment='center',horizontalalignment='center',fontsize=16)
	
	all_fluxs = []
	
	##If no repeated catalogues to combine, skip
	if num_matches==0 or num_matches==1:
		pass
	##Otherwise, plot the combined fluxes
	else:
		##Calculate and plot the combined fluxes of the two sources, even if source one or two has been accepted
		##just as a guide
		src_all = mkl.get_allinfo(all_info)
		
		if accepted_inds=='Nope':
			pass
		else:
			comb_crit, ra_ws, rerr_ws, dec_ws, derr_ws, temp_freqs, comb_freqs, comb_fluxs, comb_ferrs, comb_fit, comb_jstat, comb_chi_red, combined_names, set_freqs, set_fluxs, set_fits, comb_sources = mkl.combine_flux(src_all,src_g,accepted_inds,'plot=yes',len(matches))
		
		##If the criteria sent the double to be combined, actually plot the fitted line
		#if dom_crit == 'No dom. source':
			
			#for freq,flux in zip(set_freqs,set_fluxs):
				#ax_spectral.plot(freq,flux,linestyle='--',linewidth=1,color='r')
			
			#split_colors = ['#AE70ED','#FFB60B','#62A9FF','#59DF00']
			#for fit in set_fits:
				#ind = set_fits.index(fit)
				##ax_spectral.plot(set_freqs[ind],set_fluxs[ind],linestyle='--',linewidth=1.0,color=split_colors[ind],alpha=0.7)
				#split_p, = ax_spectral.plot(temp_freqs,np.exp(fit.params[1] + np.log(temp_freqs)*fit.params[0]),linestyle='-',linewidth=1.5,color=split_colors[ind])
				#spec_labels.append(split_p)
				#SIs.append([fit.params[0],'split %d' %(ind+1)])
			
		bright_colours = ['#FF6600','#33FF33','#FF47A3','#00ebb3']
		
		for freq in xrange(len(comb_freqs)):
			plot_errors_comb('*',bright_colours[freq],comb_freqs[freq],comb_fluxs[freq],comb_ferrs[freq],'combo',16,ax_spectral)
		comb_p, = ax_spectral.plot(temp_freqs,np.exp(comb_fit.fittedvalues),linestyle='--',linewidth=1.5,color='k')
		spec_labels.append(comb_p)
		SIs.append([comb_fit.params[0],'comb'])
		
		##Send the combined fluxes to the all_fluxs so that ax_spectral is scaled appropriately
		for flux in comb_fluxs:
			all_fluxs.append(flux)
			
			#for pos in xrange(len(ra_ws)):
				#patch = plot_pos('*',bright_colours[pos],ra_ws[pos],dec_ws[pos],rerr_ws[pos],derr_ws[pos],combined_names[pos],14,ax_main,ax_main.get_transform("fk5"))

	##==============================================================
	
	src_all2 = mkl.get_allinfo(all_info2)
		
	if accepted_inds=='Nope':
		pass
	else:
		comb_crit, ra_ws, rerr_ws, dec_ws, derr_ws, temp_freqs, comb_freqs, comb_fluxs, comb_ferrs, comb_fit, comb_jstat, comb_chi_red, combined_names, set_freqs, set_fluxs, set_fits, comb_sources = combine_flux_nobase(src_all2,src_g,xrange(len(src_all2.names)),'plot=yes',len(matches))
	
	bright_colours = ['#FF6600','#33FF33','#FF47A3','#00ebb3']
	
	for freq in xrange(len(comb_freqs)):
		plot_errors_comb('*',bright_colours[freq],comb_freqs[freq],comb_fluxs[freq],comb_ferrs[freq],'combo',16,ax_spectral_right)
	comb_p, = ax_spectral_right.plot(temp_freqs,np.exp(comb_fit.fittedvalues),linestyle='--',linewidth=1.5,color='m')
	spec_labels.append(comb_p)
	SIs.append([comb_fit.params[0],'comb_{fix}'])
	
	##Send the combined fluxes to the all_fluxs so that ax_spectral is scaled appropriately
	for flux in comb_fluxs:
		all_fluxs.append(flux)
	

	##Fill the left hand plots with information goodness
	fill_left_plots(all_info,ra_main,dec_main,ax_main,ax_spectral,tr_fk5,wcs,all_fluxs,ra_down_lim,ra_up_lim,dec_down_lim,dec_up_lim,delta_RA)
	fill_right_plots_nosource(all_info2,ra_main,dec_main,ax_main_right,ax_spectral_right,tr_fk5_right,wcs_right,all_fluxs,ra_down_lim,ra_up_lim,dec_down_lim,dec_up_lim,delta_RA)
	
	##If no repeated catalogues to combine, skip
	if num_matches==0 or num_matches==1:
		pass
	##Otherwise, plot the combined fluxes
	else:
		##Calculate and plot the combined fluxes of the two sources, even if source one or two has been accepted
		##just as a guide
		src_all = mkl.get_allinfo(all_info)
		bright_colours = ['#FF6600','#33FF33','#FF47A3','#00ebb3']
		if accepted_inds=='Nope':
			pass
		else:
			comb_crit, ra_ws, rerr_ws, dec_ws, derr_ws, temp_freqs, comb_freqs, comb_fluxs, comb_ferrs, comb_fit, comb_jstat, comb_chi_red, combined_names, set_freqs, set_fluxs, set_fits, comb_sources = mkl.combine_flux(src_all,src_g,accepted_inds,'plot=yes',len(matches))
		
			for pos in xrange(len(ra_ws)):
				patch = plot_pos_comb('*',bright_colours[pos],ra_ws[pos],dec_ws[pos],0,0,combined_names[pos],16,ax_main,ax_main.get_transform("fk5")) #rerr_ws[pos],derr_ws[pos]
				
	src_all2 = mkl.get_allinfo(all_info2)
	bright_colours = ['#FF6600','#33FF33','#FF47A3','#00ebb3']
	if accepted_inds=='Nope':
		pass
	else:
		comb_crit, ra_ws, rerr_ws, dec_ws, derr_ws, temp_freqs, comb_freqs, comb_fluxs, comb_ferrs, comb_fit, comb_jstat, comb_chi_red, combined_names, set_freqs, set_fluxs, set_fits, comb_sources = mkl.combine_flux(src_all2,src_g,xrange(len(src_all2.names)),'plot=yes',len(matches))
	
		for pos in xrange(len(ra_ws)):
			patch = plot_pos_comb('*',bright_colours[pos],ra_ws[pos],dec_ws[pos],0,0,combined_names[pos],16,ax_main_right,ax_main_right.get_transform("fk5")) #rerr_ws[pos],derr_ws[pos]
	
	#fig.tight_layout()
	fig.subplots_adjust(bottom=0.1)
	fig.subplots_adjust(left=0.15)
	
	##Make room at the top of the plot for a legend for ax_main, make the legend
	fig.subplots_adjust(top=0.85)
	
	leg_labels = [r'$\alpha_{%s}$ = %.2f' %(SI[1],SI[0]) for SI in SIs]
	main_handles,main_labels = ax_main.get_legend_handles_labels()

	main_right_handles,main_right_labels = ax_main_right.get_legend_handles_labels()
	
	for label in main_right_labels:
		if label not in main_labels and 'ombine' not in label:
			main_labels.append(label)
			main_handles.append(main_right_handles[main_right_labels.index(label)])
	
	
	main_leg = fig.add_axes([0.05,0.87,0.9,0.05])
	main_leg.axis('off')
	main_leg.legend(main_handles,main_labels,loc='lower center',prop={'size':12},ncol=8) #,bbox_to_anchor=(0,1.02),
	
	spec_leg = fig.add_axes([0.48,0.15,0.125,0.35])
	spec_leg.axis('off')

	##Stop the legend from having so many entries that it goes off the plot	
	if len(spec_labels)>11:
		trim_labels = spec_labels[:10]
		trim_labels.append(spec_labels[-1])
		trim_legs = leg_labels[:10]
		trim_legs.append(leg_labels[-1])
		spec_leg.legend(trim_labels,trim_legs,loc='center',prop={'size':14},fancybox=True)
	else:
		spec_leg.legend(spec_labels,leg_labels,loc='center',prop={'size':14},fancybox=True)
	
	##Create an axes to contain patches for an ellipse legend
	patch_leg = fig.add_axes([0.015,0.1,0.06,0.75])
	patch_leg.set_xticks([])
	patch_leg.set_yticks([])
	patch_leg.set_xticklabels([])
	patch_leg.set_yticklabels([])
	
	##See what catalogues are present in the match
	present_cats = [cat for cat in set(src_g.cats) if cat!='-100000.0']
	##Scale accordingly
	increment = 1.0/(2+len(present_cats))
	ell_positions = np.arange(increment/2,1,increment)
	
	##Find the axes coord transform
	patch_trans = patch_leg.transAxes
	##Plot and name the resolution ellipse
	ell = patches.Ellipse((0.5,ell_positions[-2]),0.9,increment-0.05,angle=0,
		transform=patch_trans, linestyle='dashed',fc='none',lw=1.1,color='gray')
	patch_leg.add_patch(ell)
	patch_leg.text(0.5,ell_positions[-2],'Resolution\n+ error',
		transform=patch_trans,verticalalignment='center',horizontalalignment='center',fontsize=14)
	
	##Plot and name the search ellipse
	ell = patches.Ellipse((0.5,ell_positions[-1]),0.9,increment-0.05,angle=0,
		transform=patch_trans, linestyle='dashdot',fc='none',lw=1.1,color='k')
	patch_leg.add_patch(ell)
	patch_leg.text(0.5,ell_positions[-1],'Search\nradius',
		transform=patch_trans,verticalalignment='center',horizontalalignment='center',fontsize=14)
	
	##Use the same method as plot_all - for some reason was getting transform errors.
	##so do it separately here (sigh)
	for cat in present_cats:
		col_ind = matched_cats.index(cat)
		position_ind = present_cats.index(cat)
		patch_leg.errorbar(0.5,ell_positions[position_ind],0.01,0.075,marker=markers[col_ind],ms=8,mfc=marker_colours[col_ind],
			mec=marker_colours[col_ind],ecolor=marker_colours[col_ind],markeredgewidth=1,label='meh',linestyle='None',transform=patch_trans)
		
		ell = patches.Ellipse((0.5,ell_positions[position_ind]),0.9,increment-0.05,angle=0, transform=patch_trans,
			fc=ell_colours1[col_ind],color=ell_colours2[col_ind],alpha=alphas[col_ind])
		patch_leg.add_patch(ell)
		
		if cat == 'nvss': cat = 'NVSS'
		elif cat == 'vlssr': cat = 'VLSSr'
		elif cat == 'comp_v11': cat = 'KGS'
		
		patch_leg.text(0.5,ell_positions[position_ind]-(increment/2-0.04),cat,
			transform=patch_trans,verticalalignment='center',horizontalalignment='center',fontsize=16)
	
	#plt.show()
	#src_all = mkl.get_allinfo(all_info)
	
	
	
	fig.savefig('%s_edit.png' %(src_all2.names[0]),bbox_inches='tight',dpi=100)
	#fig.savefig('paper_edit_top.pdf',bbox_inches='tight',dpi=200)
	
	#subprocess.call('pdftops -eps paper_edit_top.pdf paper_edit_top.eps',shell=True)
	#subprocess.call('rm paper_edit_top.pdf',shell=True)
	
	
	plt.close()

#!/usr/bin/python
import numpy as np
import statsmodels.api as sm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Ellipse
from itertools import combinations
import atpy
import optparse
import pyfits as pf
from astLib import astWCS, astPlots
import os
from matplotlib import patches
from sys import path
#path.append('/home/jline/Documents/cataloguing/PUMA/scripts')
import plot_outcomes_lib as pol
import make_table_lib as mkl
from pylab import MaxNLocator
import subprocess
import aplpy

##Set up a load of markers, colours and alpha values
markers = ['o','*','s','^','D','8','H','>','<','8','v','d']
marker_sizes = np.array([8,11,8,8,7,10,11,11,11,11,11,11,11]) + 60
marker_colours = ['m', "b", "#e5e500", "r", "g", 'k', "#c0531f",'#660066','#000099','y','#990000','#003300']
ell_colours1 = ["#C370C8", "#698ACD", "#AA6600", "#D84C77", "#5FC34F", 'k', "#D46733",'m','b','y','r','g']
ell_colours2 = ['m', "b", "#8B5A00", "r", "g", 'k', "#c0531f",'#660066','#000099','y','#990000','#003300']
alphas = [0.4,0.4,0.4,0.4,0.4,0.4,0.3,0.45,0.5,0.4,0.5,0.4]


dr= np.pi/180.0

parser = optparse.OptionParser()

parser.add_option('-f', '--cat_freqs', 
	help='Enter the frequencies of each catalogue')

parser.add_option('-m', '--matched_cats', 
	help='Enter the names of each catalogue, base catalogue first')

parser.add_option('-i', '--input_bayes', 
	help='Enter name of eyeball bayes file')

parser.add_option('-o', '--prob_thresh',
	help='The lower and upper probability thresholds - separate with a comma')

parser.add_option('-e', '--epsilon_thresh', 
	help='Cut-off threshold for the epsilon residuals')

parser.add_option('-n', '--extended_names', 
	help='List of source names to plot')

parser.add_option('-x', '--chi_thresh', 
	help='Cut-off threshold for the chi squared residuals')

parser.add_option('-r', '--resolution', 
	help='Resolution of base catalogue in "deg:arcmin:arcsec" ie "00:03:00" for 3arcmins')

parser.add_option('-b', '--mrc', default=False, action='store_true',
	help='Plot with mrc as the base rather than mwacs')

parser.add_option('-d', '--output_dir',
	help='Where to output all the plots and such')

options, args = parser.parse_args()

cat_fs = options.cat_freqs.split(',')

matched_cats = options.matched_cats.split(',')

num_freqs = []
for freq in cat_fs: num_freqs.append(len(freq.split('~')))

low_prob,high_prob = map(float,options.prob_thresh.split(','))
jstat_thresh = float(options.epsilon_thresh)
chi_thresh = float(options.chi_thresh)
closeness = mkl.dec_to_deg(options.resolution)/2

mkl.closeness = closeness
mkl.high_prob = high_prob
mkl.low_prob = low_prob
mkl.chi_thresh = chi_thresh
mkl.jstat_thresh = jstat_thresh
mkl.num_freqs = num_freqs
mkl.split = '00:02:20'

pol.closeness = closeness
pol.high_prob = high_prob
pol.low_prob = low_prob
pol.chi_thresh = chi_thresh
pol.jstat_thresh = jstat_thresh
pol.num_freqs = num_freqs
pol.split = '00:02:20'
pol.matched_cats = matched_cats



def plot_all(cat,name,ra,rerr,dec,derr,major,minor,PA,ax):
	##Plot the colour by index of catalogue in matched_cats 
	ind = matched_cats.index(cat)
	ax.show_markers(ra,dec,marker=markers[ind],s=marker_sizes[ind],facecolor=marker_colours[ind],edgecolor=marker_colours[ind],linewidth=1,label=name)
	
	#if cat=='mrc':
		#pass
	#else:
	##If no minor/major information, don't plot
	if float(minor)!=-100000.0:
		if float(major)!=-100000.0:
			ax.show_ellipses(ra,dec,float(minor),float(major),float(PA),linewidth=1.5,edgecolor=ell_colours2[ind],facecolor='none')

def find_extra(cat_name,plot_cat,cat_data,image_cats,over_plots,source_names,handles,plotted_names,ra_main,dec_main,extra_lines,freq):
	'''Finds any sources in the specified catalogue that wasn't in the original search, and plots
	them to the relevant over_plot subplot'''
	plot_ind = image_cats.index(cat_name)
	##Select the correct plot and wsc
	extra_plot = over_plots[plot_ind]
	
	##Search the catalogue and plot the found sources
	for line in cat_data[:-1]:
		info = line.split('|')[1:-1]
		#print info
		name,ra,rerr,dec,derr,flux,ferr,major,minor,PA,flag,ID = info
		name = name.split()[0]
		##Gets around long strings 
		ra,dec = float(ra[:10]),float(dec[:10])
		if (name not in source_names) and mkl.arcdist(ra_main,ra,dec_main,dec) < 4*closeness:
			flag = '-100000.0'
			ID = '-100000.0'
			
			edit_nums = [rerr,derr,flux,ferr,major,minor,PA]
			
			for i in xrange(len(edit_nums)):
				try: 
					if 'E' in edit_nums[i]:
						print edit_nums[i]
						edit_nums[i] = float(edit_nums[i][:-3]) * 10 ** float(edit_nums[i][:-2])
						print edit_nums[i]
					else:
						edit_nums[i] = float(edit_nums[i])
				except:
					edit_nums[i] = -100000.0
			
			rerr,derr,flux,ferr,major,minor,PA = edit_nums
			
			extra_line = "%s %s %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %s %s\n" %(plot_cat,name.split()[0],ra, rerr, dec, derr, freq, flux, ferr, major, major, PA,flag,ID)
			extra_lines.append(extra_line)
			
			if plot_cat == 'gleam':
				ind = 0
			else:
				ind = matched_cats.index(plot_cat)

			extra_plot.show_markers(float(ra),float(dec),marker='x',s=50,facecolor=marker_colours[ind],edgecolor=marker_colours[ind],linewidth=2,label=name)
			if float(minor)!=-100000.0:
				if float(major)!=-100000.0:
					extra_plot.show_ellipses(ra,dec,float(minor),float(major),float(PA),linewidth=1.5,edgecolor=marker_colours[ind],facecolor='none',alpha=0.6)
			
def do_plot_image(all_info,image_cats,image_files,present_cats,source_names,src_g,matches):
	##Work out how many plots across we need
	num_plots = len(image_files)
	
	fig_width = num_plots*5
	fig = plt.figure(figsize = (fig_width,10))
	
	##Find out the ra,dec of the base catalogue source
	info_base = all_info[0].split()
	name_main = info_base[1]
	ra_main = float(info_base[2])
	rerr = float(info_base[3])
	dec_main = float(info_base[4])
	derr = float(info_base[5])
	
	plot_lim = 8*closeness
	
	##Need the axes and world coordinate systems for each catalogue image to
	##plot on later
	over_plots = []
	handles = []
	plotted_names = []
	
	for cat in image_cats:
		im_ind = image_cats.index(cat)
		image_name = image_files[im_ind]
		
		ax_overplot = aplpy.FITSFigure(image_name,figure=fig, subplot=(2,num_plots,im_ind+1))
		ax_image = aplpy.FITSFigure(image_name,figure=fig, subplot=(2,num_plots,num_plots + im_ind + 1))
		
		ax_overplot.show_colorscale(cmap='gray')
		
		ax_overplot.show_ellipses(ra_main,dec_main,4*closeness,4*closeness,0,linestyle='dashdot',linewidth=1.1,edgecolor='w',facecolor='none',alpha=0.6)
		ax_overplot.show_ellipses(ra_main,dec_main,2*(closeness+rerr),2*(closeness+derr),0,linestyle='dashed',linewidth=1.1,edgecolor='w',facecolor='none',alpha=0.6)
		
		ax_image.show_colorscale(cmap='gnuplot')

		for f in [ax_overplot,ax_image]:
			f.recenter(ra_main,dec_main, width=plot_lim, height=plot_lim)
			#f.grid.set_yspacing(10.0)
			#f.add_grid()
			#f.grid.set_linewidth(2.0)
			#f.grid.set_linestyle('dashed')
			#f.grid.show()
			#f.grid.set_xspacing(0.05)  # degrees
			f.axis_labels.set_font(size=14)
			f.tick_labels.set_font(size=12)
			f.ticks.set_xspacing(0.0625)
			f.tick_labels.set_xformat('hh:mm:ss')
			f.tick_labels.set_yformat('dd:mm')
			f.add_colorbar()
			f.colorbar.set_font(size=14)
			
			f.axis_labels.hide_y()
			f.axis_labels.hide_x()
			
			if im_ind > 0: 
				f.tick_labels.hide_y()
				
			if cat == 'vlssr':
				f.add_label(0.05,0.93,'VLSSr 74MHz',color='w',weight='semibold',size='large',relative=True,horizontalalignment='left')
			elif cat == 'nvss':
				f.add_label(0.05,0.93,'NVSS 1400MHz',color='w',weight='semibold',size='large',relative=True,horizontalalignment='left')
			elif cat == 'sumss':
				f.add_label(0.05,0.93,'SUMSS 843MHz',color='w',weight='semibold',size='large',relative=True,horizontalalignment='left')
			elif cat == 'gleam':
				f.add_label(0.05,0.93,'IDR2 200MHz',color='w',weight='semibold',size='large',relative=True,horizontalalignment='left')	
				
			
		over_plots.append(ax_overplot)
	##Get the information and plot the positions on the over_plots
	#ras = []
	#decs = []
	###Use these to put in the legend later
	handles = []
	plotted_names = []
	
	for i in xrange(len(all_info)):
		info=all_info[i].split()
		cat = info[0]
		name = info[1]
		ra = float(info[2])
		#ras.append(ra)
		rerr = float(info[3])
		dec = float(info[4])
		#decs.append(dec)
		derr = float(info[5])
		major = info[-5]
		minor = info[-4]
		PA = info[-3]
		ID = info[-1]
		
		###Plot positions and elliptical fits on bottom left plot
		##plot_all(cat,name,ra,rerr,dec,derr,major,minor,PA,ax_main)
		###If catalogue in one of the images, plot the informaion on the applicaable over_plot
				
		if options.mrc:
			if cat=='mrc':
				for i in xrange(len(over_plots)): 
					plot_all(cat,name,ra,rerr,dec,derr,major,minor,PA,over_plots[i])
				
		else:
			if 'gleam' in cat:
				#pass
				for i in xrange(len(over_plots)): 
					plot_all(cat,name,ra,rerr,dec,derr,major,minor,PA,over_plots[i])
				
		if cat in image_cats:
			cat_ind = image_cats.index(cat)
			plot_all(cat,name,ra,rerr,dec,derr,major,minor,PA,over_plots[cat_ind])
			
	extra_lines = []

	###Plot cat entries that weren't in the original match
	if 'sumss' in image_cats:
		find_extra('sumss','gleam',gleam_data,image_cats,over_plots,source_names,handles,plotted_names,ra_main,dec_main,extra_lines,200.0)
		find_extra('sumss','sumss',sumss_data,image_cats,over_plots,source_names,handles,plotted_names,ra_main,dec_main,extra_lines,843.0)
		
	if 'vlssr' in image_cats:
		find_extra('vlssr','gleam',gleam_data,image_cats,over_plots,source_names,handles,plotted_names,ra_main,dec_main,extra_lines,200.0)
		find_extra('vlssr','vlssr',vlssr_data,image_cats,over_plots,source_names,handles,plotted_names,ra_main,dec_main,extra_lines,74.0)
		
	if 'nvss' in image_cats:
		find_extra('nvss','gleam',gleam_data,image_cats,over_plots,source_names,handles,plotted_names,ra_main,dec_main,extra_lines,200.0)
		find_extra('nvss','nvss',nvss_data,image_cats,over_plots,source_names,handles,plotted_names,ra_main,dec_main,extra_lines,1400.0)
		
	if 'gleam' in image_cats:
		find_extra('gleam','gleam',gleam_data,image_cats,over_plots,source_names,handles,plotted_names,ra_main,dec_main,extra_lines,200.0)
		
	extras_file = open("%s/%s_extra-sources.txt" %(options.output_dir,source_names[0]) ,'w+')

	if len(extra_lines) > 0:
		for line in extra_lines:
			extras_file.write(line)

	
	###Make room at the top of the plot for a legend for ax_main, make the legend
	fig.subplots_adjust(top=0.85)
	
	for over_plot in over_plots:
		ax_handles,ax_labels = over_plot._ax1.get_legend_handles_labels()
		for handle,label in zip(ax_handles,ax_labels):
			handles.append(handle)
			plotted_names.append(label)
	
	leg_handles = []
	leg_labels = []
	
	for hand,lab in zip(handles,plotted_names):
		if lab not in leg_labels:
			leg_handles.append(hand)
			leg_labels.append(lab)
			
	main_leg = fig.add_axes([0.05,0.87,0.9,0.05])
	main_leg.axis('off')

	main_leg.legend(leg_handles,leg_labels,loc='lower center',prop={'size':11},ncol=2*num_plots,scatterpoints=1) #,bbox_to_anchor=(0,1.02),
	
	return fig
		
##======================================================================================
##Check whether the source is already a HACK

##Catalogue data
simple_cats = "./simple_cats"
mrc_data = open('%s/mrc_simple.txt' %simple_cats).read().split('\n')
sumss_data = open('%s/sumss_simple.txt' %simple_cats).read().split('\n')
vlssr_data = open('%s/vlssr_simple.txt' %simple_cats).read().split('\n')
nvss_data = open('%s/nvss_simple.txt' %simple_cats).read().split('\n')
mwacs_data = open('%s/mwacs_simple.txt' %simple_cats).read().split('\n')
gleam_data = open('%s/IDR2_noneg_good.txt' %simple_cats).read().split('\n')
	
##======================================================================================

##Input eyeball file
bayes_comp = open(options.input_bayes).read().split('END_GROUP')
del bayes_comp[-1]

extended_lines = open(options.extended_names).read().split('\n')

if extended_lines[-1] == '': del extended_lines[-1]

extended_names = [line.split()[0] for line in extended_lines]

fail = open('%s/failed_runs.txt' %options.output_dir,'w+')

#print len(bayes_comp)

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
	src_g = mkl.get_srcg(match1)
	##Get some info and find which catalogues are present 
	src_all = mkl.get_allinfo(all_info)
	
	try:
		if src_all.names[0] in extended_names:
			present_cats = [cat for cat in src_all.cats if cat!='-100000.0']
			
			meh,num_matches,accept_matches,accepted_inds,accept_type,stage = stats.split()
			
			image_num = bayes_comp.index(comp)+1
			image_locs = "./extended_fits"
			image_cats = [cat for cat in ['gleam','vlssr','sumss','nvss'] if os.path.exists("%s/%s_%s.fits" %(image_locs,src_all.names[0],cat))==True]
			image_files =["%s/%s_%s.fits" %(image_locs,src_all.names[0],cat) for cat in ['gleam','vlssr','sumss','nvss'] if os.path.exists("%s/%s_%s.fits" %(image_locs,src_all.names[0],cat))==True]
			
			fig = do_plot_image(all_info,image_cats,image_files,present_cats,src_all.names,src_g,matches)
			#plt.tight_layout()
			fig.savefig('%s/%s_extended.png' %(options.output_dir,src_all.names[0]),bbox_inches='tight',dpi=100)
			plt.close()
			
	except:
		print src_all.names[0], ' failed'
		fail.write(src_all.names[0]+'\n')
		
fail.close()

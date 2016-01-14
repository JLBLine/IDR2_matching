#!/usr/bin/python
import atpy
import numpy as np
import plot_outcomes_lib_mod as pol
from sys import path
#path.append('/home/jline/Documents/cataloguing/FHD/carroll/paper/complete_v11/investigate_possibles_fix')
import make_table_lib_mod as mkl
import optparse
import copy

parser = optparse.OptionParser()

parser.add_option('-p', '--pref_cats',
	help='Enter names of matched cataloges, in order of preferable position')
	
parser.add_option('-m', '--matched_cats',
	help='Enter names of matched cataloges, in order of match')
	
parser.add_option('-c', '--cat_freqs', 
	help='Enter number of frequencies in each catalogue')
	
parser.add_option('-i', '--input_bayes', 
	help='Enter name of input matched bayes')

parser.add_option('-f', '--input_fiddle', 
	help='Enter name of input modified bayes')

parser.add_option('-o', '--prob_thresh',
	help='The lower and upper probability thresholds - separate with a comma')

parser.add_option('-e', '--epsilon_thresh', 
	help='Cut-off threshold for the epsilon residuals')

parser.add_option('-x', '--chi_thresh', 
	help='Cut-off threshold for the chi squared residuals')

parser.add_option('-r', '--resolution', 
	help='Resolution of base catalogue in "deg:arcmin:arcsec" ie "00:03:00" for 3arcmins')

parser.add_option('-s', '--split',default=0, 
	help='The resolution ("deg:arcmin:arcsec") over which to split combined sources')

parser.add_option('-n', '--choices',
	help='Names and choices of sources to plot')

parser.add_option('-t', '--tag',
	help='Names tag of the plot')

parser.add_option('-z', '--combine_all_only',default=False,action='store_true',
	help='Enable if plotting the combine all "c" option used by use_choices.py')

options, args = parser.parse_args()

##Set up a bunch of initial parameters-----------------------------------------
##-----------------------------------------------------------------------------
low_prob,high_prob = map(float,options.prob_thresh.split(','))
jstat_thresh = float(options.epsilon_thresh)
chi_thresh = float(options.chi_thresh)
closeness = mkl.dec_to_deg(options.resolution)/2

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
split = options.split
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------

##Transfer variables to the mkl and pol modules
pol.closeness = closeness
pol.high_prob = high_prob
pol.low_prob = low_prob
pol.chi_thresh = chi_thresh
pol.jstat_thresh = jstat_thresh
pol.num_freqs = num_freqs
pol.split = split
pol.matched_cats = matched_cats

mkl.closeness = closeness
mkl.high_prob = high_prob
mkl.low_prob = low_prob
mkl.chi_thresh = chi_thresh
mkl.jstat_thresh = jstat_thresh
mkl.num_freqs = num_freqs
mkl.split = split

#names = open(options.names,'r').read().split('\n')
#if names[-1]=='': del names[-1]

#names = [entry.split()[0] for entry in names if entry.split()[1]=='r']


##Open the input text file (output from calculate_bayes.py)
bayes_comp = open(options.input_bayes).read().split('END_GROUP')
del bayes_comp[-1]

from matplotlib import rc
font = {'size': 14}
rc('font', **font)

choices = open(options.choices,'r').read().split('\n')
if choices[-1] == '': del choices[-1]

choice_dict = {}
choice_names = []

if options.combine_all_only:
	for choice in choices:
		try:
			name,inc_names = choice.split()
			#if inc_names == 'c':
				#choice_dict[name] = ['all']
				#choice_names.append(name)
			#else:
			if inc_names != 'r':
				choice_dict[name] = inc_names
				choice_names.append(name)
		##This means there are 3 input, ie name a c
		except:
			pass

else:
	bayes_fiddle = open(options.input_fiddle).read().split('END_GROUP')
	if bayes_fiddle[-1]=='': del bayes_fiddle[-1]
	if bayes_fiddle[-1].split('\n')==['','','']: del bayes_fiddle[-1]
	if bayes_fiddle[-1].split('\n')==['','']: del bayes_fiddle[-1]
	for choice in choices:
		name,inc_names = choice.split()
		inc_names = inc_names.split(',')
		choice_dict[name] = inc_names
		choice_names.append(name)

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
	
	meh,num_matches,accept_matches,accepted_inds,accept_type,stage = stats.split()

	match_crit = 'meh'
	##This line applies positional criteria, and tells us if a simple one catalogue repeat source, returning
	##how many combinations are possible and statistics on them
	match_crit = "%d of %d \ncombinations \npossible" %(int(accept_matches),len(matches))
	
	
	##Put the information for every source in the matched group in one source_group() class
	##(see apply_criteria_lib for source_group())
	src_all = mkl.get_allinfo(all_info)
	
	if src_all.names[0] in choice_names: #and src_all.names[0]=='8514.0':
		##If just combining all infomation present, use the same bayes file
		if options.combine_all_only:
			
			chunks2 = comp.split('START_COMP')
			
			all_info2 = chunks2[0].split('\n')
			
			for entry in all_info2:   
				if entry=='': del all_info2[all_info2.index(entry)]
			for entry in all_info2:
				if 'START' in entry: del all_info2[all_info2.index(entry)]
			
			src_all2 = mkl.get_allinfo(all_info2)
			
			if type(accepted_inds)==list:
					pass
			else:
				accepted_inds = map(int,accepted_inds.split(','))
				
			user_edit = choice_dict[src_all.names[0]]
				
			if user_edit == 'c':
				user_edit = ['all']
				#accept_matches2 = accept_matches
				#accepted_inds2 = accepted_inds
				this_match = 'Nope'
			else:
				accepted_match = matches[int(user_edit)-1].split()
				this_match = int(user_edit)
				user_edit = [name for name in mkl.get_srcg(accepted_match).names if name != '-100000.0']
				
				#repeats = [src_all.names[i] for i in xrange(len(src_all.names)) if src_all.cats.count(src_all.cats[i]) > 1]
				##these_names = 
				#accepted_inds2 = [src_all.names.index(name) for name in repeats]
				#accept_matches2 = 'waaah'
				#accept_matches2 = accept_matches
				
			#print accepted_inds,accepted_inds2


			pol.create_plot_nocombo(comp,comp,accepted_inds,int(accept_matches),this_match,match_crit,stage,accept_type,options.tag,user_edit)
		##Otherwise, need to get the fiddled bayes info
		else:
			for comp2 in bayes_fiddle:
				chunks2 = comp2.split('START_COMP')
				
				all_info2 = chunks2[0].split('\n')
				
				for entry in all_info2:   
					if entry=='': del all_info2[all_info2.index(entry)]
				for entry in all_info2:
					if 'START' in entry: del all_info2[all_info2.index(entry)]
				
				src_all2 = mkl.get_allinfo(all_info2)
					
				if src_all2.names[0] == src_all.names[0]:

					if type(accepted_inds)==list:
						pass
					else:
						accepted_inds = map(int,accepted_inds.split(','))
					
					pol.create_plot_nocombo(comp,comp2,accepted_inds,int(accept_matches),match_crit,stage,accept_type,options.tag,choice_dict[src_all.names[0]])
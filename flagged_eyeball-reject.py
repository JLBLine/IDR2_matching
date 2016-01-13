from sys import path
path.append('./PLOTS/plotting_scripts')
import make_table_lib as mkl
#import optparse

reject_file = open('IDR2_reject_sources.txt','w+')
eyeball_file = open('IDR2_eyeball_sources.txt','w+')

reject_list = []
eyeball_list = []

bayes_comp = open('./multifreq/puma_gleammulti-v-m-s-n-a-eyeball.txt').read().split('END_GROUP')
del bayes_comp[-1]

for comp in bayes_comp:

	##Get the information into nice usable forms, and get rid of empty/pointless
	##entries
	chunks = comp.split('START_COMP')
	all_info = chunks[0].split('\n')
	
	for entry in all_info:   
		if entry=='': del all_info[all_info.index(entry)]
	for entry in all_info:
		if 'START' in entry: del all_info[all_info.index(entry)]
		
	matches = chunks[1].split('\n')
	
	del matches[0],matches[-1]
	stats = matches[-1]
	del matches[-2:]
	
	##Get some info and find which catalogues are present 
	src_all = mkl.get_allinfo(all_info)

	meh,num_matches,accept_matches,accepted_inds,accept_type,stage = stats.split()

	if accept_type == 'reject':
		reject_list.append(src_all.names[0])
	else:
		eyeball_list.append(src_all.names[0])
		
bayes_comp = open('./deep/puma_gleamdeep-v-m-s-n-a-eyeball.txt').read().split('END_GROUP')
del bayes_comp[-1]

for comp in bayes_comp:

	##Get the information into nice usable forms, and get rid of empty/pointless
	##entries
	chunks = comp.split('START_COMP')
	all_info = chunks[0].split('\n')
	
	for entry in all_info:   
		if entry=='': del all_info[all_info.index(entry)]
	for entry in all_info:
		if 'START' in entry: del all_info[all_info.index(entry)]
		
	matches = chunks[1].split('\n')
	
	del matches[0],matches[-1]
	stats = matches[-1]
	del matches[-2:]
	
	##Get some info and find which catalogues are present 
	src_all = mkl.get_allinfo(all_info)

	meh,num_matches,accept_matches,accepted_inds,accept_type,stage = stats.split()

	if accept_type == 'reject':
		#print 'here'
		#if src_all.names[0] not in reject_list:
		reject_list.append(src_all.names[0])
	else:
		#if src_all.names[0] not in eyeball_list:
		eyeball_list.append(src_all.names[0])
		
for name in sorted(set(reject_list)): reject_file.write(name+'\n')
for name in sorted(set(eyeball_list)): eyeball_file.write(name+'\n')

reject_file.close()
eyeball_file.close()
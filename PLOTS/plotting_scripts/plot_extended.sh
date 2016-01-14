#!/bin/sh
# python plot_extended.py	--matched_cats=gleam_multi,vlssr,mrc,sumss,nvss,atg20 \
# 	--input_bayes=../../multifreq/puma_gleammulti-v-m-s-n-a-eyeball.txt \
# 	--cat_freqs=76~84~92~99~107~115~122~130~143~151~158~166~174~181~189~197~204~212~220~227,74,408,843,1400,5000~8000~20000 \
# 	--prob_thresh=0.8,0.95 \
# 	--epsilon_thresh=0.1 --chi_thresh=10 \
# 	--resolution=00:02:20 \
# 	--extended_names=../../IDR2_eyeball_sources.txt \
# 	--output_dir=/media/jline/TOSHIBA_jline/Brian_ftw/PLOTS/eyeball

python plot_extended.py	--matched_cats=gleam_multi,vlssr,mrc,sumss,nvss,atg20 \
	--input_bayes=../../multifreq/puma_gleammulti-v-m-s-n-a-accept.txt \
	--cat_freqs=76~84~92~99~107~115~122~130~143~151~158~166~174~181~189~197~204~212~220~227,74,408,843,1400,5000~8000~20000 \
	--prob_thresh=0.8,0.95 \
	--epsilon_thresh=0.1 --chi_thresh=10 \
	--resolution=00:02:20 \
	--extended_names=../../IDR2_eyeball_sources.txt \
	--output_dir=/media/jline/TOSHIBA_jline/IDR2_matching/PLOTS/eyeball
	
	
# 	--extended_names=../../multifreq/gleammulti_eyeball_sources.txt \
	
	

	
	
# 	--query=502.0,2271.0,2556.0,262.0,1823.0,1378.0,80.0,2449.0,2411.0,8976.0,8337.0,8063.0

#!/bin/sh
# python plot_edits.py --matched_cats=gleam_multi,vlssr,mrc,sumss,nvss,atg20 \
# 	--pref_cats=nvss,sumss,gleam_multi,atg20,mrc,vlssr \
# 	--input_bayes=../../multifreq/eyeball_reject_pile.txt \
# 	--input_fiddle=../../multifreq/eyeball_reject_pile_mod.txt \
# 	--cat_freqs=76~84~92~99~107~115~122~130~143~151~158~166~174~181~189~197~204~212~220~227,74,408,843,1400,5000~8000~20000 \
# 	--prob_thresh=0.8,0.95 \
# 	--epsilon_thresh=0.1 --chi_thresh=10 \
# 	--resolution=00:02:18 \
# 	--choices=../../multifreq/eyeball_reject_choices.txt --tag=eyeball
python plot_edits.py --matched_cats=gleam_multi,vlssr,mrc,sumss,nvss,atg20 \
	--pref_cats=nvss,sumss,gleam_multi,atg20,mrc,vlssr \
	--input_bayes=../../multifreq/puma_gleammulti-v-m-s-n-a-eyeball.txt \
	--cat_freqs=76~84~92~99~107~115~122~130~143~151~158~166~174~181~189~197~204~212~220~227,74,408,843,1400,5000~8000~20000 \
	--prob_thresh=0.8,0.95 \
	--epsilon_thresh=0.1 --chi_thresh=10 \
	--resolution=00:02:18 \
	--choices=../../multifreq/test_choices.txt --tag=eyeball --combine_all_only
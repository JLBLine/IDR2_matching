#!/bin/sh
 plot_outcomes.py --matched_cats=comp_v11,vlssr,mrc,sumss,nvss \
 	--pref_cats=nvss,sumss,comp_v11,mrc,vlssr \
 	--input_bayes=../bayes_v11+v10extras-v-m-s-n.txt \
 	--cat_freqs=182.435,74,408,843,1400 \
 	--prob_thresh=0.8,0.95 \
 	--epsilon_thresh=0.1 --chi_thresh=10 \
 	--resolution=00:02:18 \
 	--query=298.0 \
  	--write=all

# python do_plot_outcomes.py ../eyeball_sources.txt
# mv *png ./eyeball
# 
#python do_plot_outcomes.py ../multiple_sources.txt
#mv *png ./multiple
## 
#python do_plot_outcomes.py ../extreme_names_v10+v11.txt
#mv *png ./extreme

#python do_plot_outcomes.py ../reject_sources.txt
#mv *png ./reject

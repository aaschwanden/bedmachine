#!/bin/bash

set -x

R=l
GS=250
fill_value=-2000000000

alpha=0.0


for gamma in 0.0 0.1 0.5 1.0 10.0 100.0
do
    for project in "jakobshavn"
    do
        for exp in  "nodhdt_nobmelt" "dhdt_bmelt" "nodhdt_bmelt"
        do
	    for scale in 0.8 0.9 1.0
	    do
		nc2cdo.py ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc
		cdo ifnotthen ${project}_mask_${GS}m.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc ${project}/tmp_${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc
		mv ${project}/tmp_${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc
		ncks -A -v x,y,mapping ${project}_mask_${GS}m.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc
		ncatted -a grid_mapping,thk,o,c,"mapping" ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc
		extract-profile.py ~/data/data_sets/Jakobshavn1985/jakobshavn_xprofile_str_point.shp ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc ${project}/profile_str_${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc
		extract-profile.py ~/data/data_sets/Jakobshavn1985/jakobshavn_xprofile_s8_point.shp ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc ${project}/profile_s8_${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc
		extract-profile.py ~/data/data_sets/GreenlandFlightlines/jakobshavn_xprofile_20090406_points.shp ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc ${project}/profile_ctr_${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc
		extract-profile.py ~/data/data_sets/GreenlandFlightlines/jakobshavn_xprofile_20090402-02-08.shp ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc ${project}/profile_ptr_${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc
	    done
	done
    done
done


for project in "jakobshavn"
do
    for dataset in "cresis"
    do
        # thk of input data sets
        #im-plot.py -p medium --colorbar_label --bounds 0 2500 -v thk --singlerow --colormap /Users/andy/base/PyPISMTools/colormaps/sst.cpt -o ${project}/${project}_${dataset}_${GS}m_thk.pdf ${project}_${dataset}_${GS}m.nc
        # topg of input data sets
        #im-plot.py -p medium --colorbar_label -v topg --singlerow --colormap /Users/andy/base/pypismtools/colormaps/wiki-2.0.cpt -o ${project}/${project}_${dataset}_${GS}m_topg.pdf ${project}_${dataset}_${GS}m.nc
        extract-profile.py ~/data/data_sets/Jakobshavn1985/jakobshavn_xprofile_str_point.shp ${project}_${dataset}_${GS}m.nc profile_str_${project}_${dataset}_${GS}m.nc
        extract-profile.py ~/data/data_sets/Jakobshavn1985/jakobshavn_xprofile_s8_point.shp ${project}_${dataset}_${GS}m.nc profile_s8_${project}_${dataset}_${GS}m.nc
        extract-profile.py  ~/data/data_sets/GreenlandFlightlines/jakobshavn_xprofile_20090406_points.shp ${project}_${dataset}_${GS}m.nc profile_ctr_${project}_${dataset}_${GS}m.nc
        extract-profile.py  ~/data/data_sets/GreenlandFlightlines/jakobshavn_xprofile_20090402-02-08.shp ${project}_${dataset}_${GS}m.nc profile_ctr_${project}_${dataset}_${GS}m.nc
    done
done



exit

for gamma in 0.0 0.5 1.0 2.0 5.0 10.0 20.0 50.0 100.0
do
    for project in "jakobshavn"
    do
        for exp in  "nodhdt_nobmelt" "dhdt_bmelt"
        do
	    for scale in 0.7 0.8 0.9 1.0
	    do
                # bed elevation
		im-plot.py -p medium --colorbar_label -v topg --singlerow --colormap ~/base/pypismtools/colormaps/wiki-2.0.cpt -o ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_${exp}_topg.pdf ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc
               # thickness  
		im-plot.py -p medium --colorbar_label --bounds 0 2500 -v thk --singlerow --colormap ~/base/pypismtools/colormaps/sst.cpt -o ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}_thk.pdf ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc
            # thk difference
		im-plot.py -p medium --colorbar_label -v thk --singlerow --colormap RdBu_r  --bounds -100 100 --obs_file ${project}_flightlines_${GS}m.nc -o ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}_thkdiff.pdf ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc
            # divHU
		im-plot.py -p medium --colorbar_label -v divHU --singlerow --colormap RdBu_r  --bounds -100 100 -o ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}_divHU.pdf ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc
            # thk difference exp-obs
		dataset=cresis
		im-plot.py -p medium --colorbar_label -v thk --singlerow --bounds -200 200 --colormap RdBu_r -o ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}_${dataset}_thk_diff.pdf --obs_file ${project}_${dataset}_${GS}m.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc
		extract-profile.py jakobshavn_xprofile_str_point.shp ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc ${project}/profile_str_${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc
		dataset=flightlines
		im-plot.py -p medium --colorbar_label -v thk --singlerow --bounds -200 200 --colormap RdBu_r -o ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}_${dataset}_thk_diff.pdf --obs_file ${project}_${dataset}_${GS}m.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc
	    done
	done
    done
done


exit


alpha=0.0
exp=nodhdt_nobmelt
project=jakobshavn
gamma=5.0
GS=250

im-plot.py -p twocol --inner_titles 'u(1985),u(2008)' --colorbar_label -v divHU_cresis --singlerow --colormap RdBu_r  --bounds -500 500 -o ${project}_${GS}m_cresis_vscale_${scale}_${exp}_divHU.pdf ../bedmachine_250m_1985/${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc 

 im-plot.py -p medium --colorbar_label -v divHU_cresis --singlerow --colormap RdBu_r  --bounds -500 500 -o ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}_divHU_cresis.pdf ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc


alpha=0.0
exp=nodhdt_nobmelt
project=jakobshavn
gamma=5.0
GS=250


extract-profile.py jakobshavn_xprofile_str_point.shp ${project}_${dataset}_${GS}m.nc profile_str_${project}_${dataset}_${GS}m.nc


~/base/pypismtools/scripts/profile-plot.py -o profile_s8_gammas_beta_0.9 --split 2 -v topg --labels 'reference,0,0.1,0.5,1,10,100' --figure_title S8 profile_s8_jakobshavn_cresis_250m.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_0.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_0.1_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/profile_s8_jakobshavn_cresis_250m.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_0.0_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_0.1_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.9_nodhdt_bmelt.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_str_gammas_beta_0.9 --split 2 -v topg --labels 'reference,0,0.1,0.5,1,10,100' --figure_title STR profile_str_jakobshavn_cresis_250m.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_0.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_0.1_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/profile_str_jakobshavn_cresis_250m.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_0.0_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_0.1_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.9_nodhdt_bmelt.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_ctr_gammas_beta_0.9 --split 2 -v topg --labels 'reference,0,0.1,0.5,1,10,100' --figure_title CTR profile_ctr_jakobshavn_cresis_250m.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_0.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_0.1_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/profile_ctr_jakobshavn_cresis_250m.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_0.0_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_0.1_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.9_nodhdt_bmelt.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_ptr_gammas_beta_0.9 --split 2 -v topg --labels 'reference,0,0.1,0.5,1,10,100' --figure_title PTR profile_ptr_jakobshavn_cresis_250m.nc jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_0.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_0.1_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/profile_ptr_jakobshavn_cresis_250m.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_0.0_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_0.1_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.9_nodhdt_bmelt.nc


~/base/pypismtools/scripts/profile-plot.py -o profile_s8_betas_gamma_0.5 --split 2 -v topg --labels 'reference,0.8,0.9,1.0' --figure_title S8 profile_s8_jakobshavn_cresis_250m.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.8_nodhdt_bmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_1.0_nodhdt_bmelt.nc profile_s8_jakobshavn_cresis_250m.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.8_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_1.0_nodhdt_bmelt.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_str_betas_gamma_0.5 --split 2 -v topg --labels 'reference,0.8,0.9,1.0' --figure_title STR profile_str_jakobshavn_cresis_250m.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.8_nodhdt_bmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_1.0_nodhdt_bmelt.nc profile_str_jakobshavn_cresis_250m.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.8_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_1.0_nodhdt_bmelt.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_ctr_betas_gamma_0.5 --split 2 -v topg --labels 'reference,0.8,0.9,1.0' --figure_title CTR profile_ctr_jakobshavn_cresis_250m.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.8_nodhdt_bmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_1.0_nodhdt_bmelt.nc profile_ctr_jakobshavn_cresis_250m.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.8_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_1.0_nodhdt_bmelt.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_ptr_betas_gamma_0.5 --split 2 -v topg --labels 'reference,0.8,0.9,1.0' --figure_title PTR profile_ptr_jakobshavn_cresis_250m.nc jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.8_nodhdt_bmelt.nc jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_1.0_nodhdt_bmelt.nc profile_ptr_jakobshavn_cresis_250m.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.8_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_1.0_nodhdt_bmelt.nc


~/base/pypismtools/scripts/profile-plot.py -o profile_s8_betas_gamma_1.0 --split 2 -v topg --labels 'reference,0.8,0.9,1.0' --figure_title S8 profile_s8_jakobshavn_cresis_250m.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_bmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_1.0_nodhdt_bmelt.nc profile_s8_jakobshavn_cresis_250m.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_1.0_nodhdt_bmelt.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_str_betas_gamma_1.0 --split 2 -v topg --labels 'reference,0.8,0.9,1.0' --figure_title STR profile_str_jakobshavn_cresis_250m.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_bmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_1.0_nodhdt_bmelt.nc profile_str_jakobshavn_cresis_250m.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_1.0_nodhdt_bmelt.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_ctr_betas_gamma_1.0 --split 2 -v topg --labels 'reference,0.8,0.9,1.0' --figure_title CTR profile_ctr_jakobshavn_cresis_250m.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_bmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_1.0_nodhdt_bmelt.nc profile_ctr_jakobshavn_cresis_250m.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_1.0_nodhdt_bmelt.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_ptr_betas_gamma_1.0 --split 2 -v topg --labels 'reference,0.8,0.9,1.0' --figure_title PTR profile_ptr_jakobshavn_cresis_250m.nc jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_bmelt.nc jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_1.0_nodhdt_bmelt.nc profile_ptr_jakobshavn_cresis_250m.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_nodhdt_bmelt.nc ../bedmachine_250m_1985_cresis/jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_1.0_nodhdt_bmelt.nc




project=jakobshavn
GS=250
alpha=0.0
exp=nodhdt_bmelt
for scale in 0.9;
do
    im-plot.py -p small_font --colorbar_position right --inner_titles '$\gamma$=0,$\gamma$=0.1,$\gamma$=0.5,$\gamma$=1,$\gamma$=10,$\gamma$=100' --colorbar_label -v topg --colormap /Users/andy/base/pypismtools/colormaps/wiki-2.0.cpt -o ${project}_${GS}m_alpha_${alpha}_gammas_vscale_${scale}_${exp}_topg.pdf ${project}/${project}_${GS}m_alpha_${alpha}_gamma_0.0_vscale_${scale}_${exp}.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_0.1_vscale_${scale}_${exp}.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_0.5_vscale_${scale}_${exp}.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_1.0_vscale_${scale}_${exp}.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_10.0_vscale_${scale}_${exp}.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_100.0_vscale_${scale}_${exp}.nc
    im-plot.py -p small_font --colorbar_position right --inner_titles '$\gamma$=0,$\gamma$=0.1,$\gamma$=0.5,$\gamma$=1,$\gamma$=10,$\gamma$=100' --colorbar_label -v divHU --colormap RdBu_r  --bounds -500 500 -o ${project}_${GS}m_alpha_${alpha}_gammas_vscale_${scale}_${exp}_divHU.pdf ${project}/${project}_${GS}m_alpha_${alpha}_gamma_0.0_vscale_${scale}_${exp}.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_0.1_vscale_${scale}_${exp}.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_0.5_vscale_${scale}_${exp}.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_1.0_vscale_${scale}_${exp}.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_10.0_vscale_${scale}_${exp}.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_100.0_vscale_${scale}_${exp}.nc
    im-plot.py -p small_font --colorbar_position right --inner_titles '$\gamma$=0,$\gamma$=0.1,$\gamma$=0.5,$\gamma$=1,$\gamma$=10,$\gamma$=100' --colorbar_label -v res_flux --colormap RdBu_r  --bounds -500 500 -o ${project}_${GS}m_alpha_${alpha}_gammas_vscale_${scale}_${exp}_res_flux.pdf ${project}/${project}_${GS}m_alpha_${alpha}_gamma_0.0_vscale_${scale}_${exp}.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_0.1_vscale_${scale}_${exp}.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_0.5_vscale_${scale}_${exp}.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_1.0_vscale_${scale}_${exp}.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_10.0_vscale_${scale}_${exp}.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_100.0_vscale_${scale}_${exp}.nc
done

~/base/pypismtools/scripts/profile-plot.py -o profile_s8_dhdt_beta_0.9 --split 2 -v topg --labels 'reference,0,0.1,0.5,1,10,100' --figure_title S8 profile_s8_jakobshavn_cresis_250m.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_0.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_0.1_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.9_nodhdt_bmelt.nc profile_s8_jakobshavn_cresis_250m.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_0.0_vscale_0.9_dhdt_bmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_0.1_vscale_0.9_dhdt_bmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_dhdt_bmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_dhdt_bmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.9_dhdt_bmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.9_dhdt_bmelt.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_str_dhdt_beta_0.9 --split 2 -v topg --labels 'reference,0,0.1,0.5,1,10,100' --figure_title STR profile_str_jakobshavn_cresis_250m.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_0.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_0.1_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.9_nodhdt_bmelt.nc profile_str_jakobshavn_cresis_250m.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_0.0_vscale_0.9_dhdt_bmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_0.1_vscale_0.9_dhdt_bmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_dhdt_bmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_dhdt_bmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.9_dhdt_bmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.9_dhdt_bmelt.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_ctr_dhdt_beta_0.9 --split 2 -v topg --labels 'reference,0,0.1,0.5,1,10,100' --figure_title CTR profile_ctr_jakobshavn_cresis_250m.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_0.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_0.1_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.9_nodhdt_bmelt.nc profile_ctr_jakobshavn_cresis_250m.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_0.0_vscale_0.9_dhdt_bmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_0.1_vscale_0.9_dhdt_bmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_dhdt_bmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_dhdt_bmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.9_dhdt_bmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.9_dhdt_bmelt.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_ptr_dhdt_beta_0.9 --split 2 -v topg --labels 'reference,0,0.1,0.5,1,10,100' --figure_title PTR profile_ptr_jakobshavn_cresis_250m.nc jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_0.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_0.1_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.9_nodhdt_bmelt.nc jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.9_nodhdt_bmelt.nc profile_ptr_jakobshavn_cresis_250m.nc jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_0.0_vscale_0.9_dhdt_bmelt.nc jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_0.1_vscale_0.9_dhdt_bmelt.nc jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_dhdt_bmelt.nc jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_dhdt_bmelt.nc jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.9_dhdt_bmelt.nc jakobshavn/profile_ptr_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.9_dhdt_bmelt.nc

im-plot.py -p twocol --colorbar_label -v dhdt  --bounds -100 100 --colormap RdBu_r --singlerow -o dhdt.pdf jakobshavn_dhdt_250m.nc

im-plot.py -p twocol --colorbar_label -v topg  --colormap /Users/andy/base/pypismtools/colormaps/wiki-2.0.cpt --singlecolumn -o dhdt_comparison_topg.pdf jakobshavn/jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_nodhdt_bmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_dhdt_bmelt.nc

im-plot.py -p twocol --colorbar_label -v divHU  --bounds -500 500 --colormap RdBu_r --singlecolumn -o dhdt_comparison_divHU.pdf jakobshavn/jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_nodhdt_bmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_dhdt_bmelt.nc



exit













im-plot.py -p twocol --colorbar_label -v magnitude --inner_titles '1985,2008' --singlerow --colormap /Users/andy/base/pypismtools/colormaps/Full_saturation_spectrum_CCW_orange.cpt -o speed_1985_2008.pdf jakobshavn_surf_vels_1985_250m.nc jakobshavn_surf_vels_2008_2009_250m.nc

im-plot.py --bounds -1000 2100 -p twocol --colorbar_label -v usurf --colormap ~/base/pypismtools/colormaps/wiki-2.0.cpt --inner_titles '1985,2008' --singlerow -o usurf_1985_2008.pdf jakobshavn_usurf_1985_250m.nc jakobshavn_usurf_2008_250m.nc

im-plot.py -p twocol --colorbar_label -v topg --inner_titles 'Bamber 2001,CReSIS' --singlerow --colormap ~/base/pypismtools/colormaps/wiki-2.0.cpt -o topg_bamber_cresis_topg_250m.pdf jakobshavn_searise_v1.1_250m.nc jakobshavn_cresis_250m.nc

gdalwarp -overwrite -r bilinear -te -220000 -2300000 -120000 -2240000 -tr 100 100 -t_srs EPSG:3413 /Volumes/data/data_sets/GreenlandLandsat/LC80090112013101LGN01_pansharpended_natural_z.tif kansas/jakobshavn_landsat8.tif

basemap-plot.py --geotiff_file jakobshavn_gimpdem_90m_hillshade.tif  -p twocol --colorbar_label -v topg --inner_title --singlerow --colormap ~/base/pypismtools/colormaps/wiki-2.0.cpt -o topg_bamber_cresis_topg_250m.pdf Greenland_5km_v1.1.nc jakobshavn_cresis_250m.nc

~/base/pypismtools/scripts/basemap-plot.py --geotiff_file jakobshavn_gimpdem_90m_hillshade.tif --alpha 0.5 -p twocol --colorbar_label -v topg --inner_titles 'Bamber (2001),Bamber (2013),CReSIS (2012)' --singlecol --colormap ~/base/pypismtools/colormaps/wiki-2.0.cpt -o bamber_cresis_topg_250m.pdf Greenland_5km_v1.1.nc ~/data/data_sets/greenland_beds_v2/Greenland_1km_v2.nc jakobshavn_cresis_250m.nc


extract-profile.py ~/data/data_sets/GreenlandFlightlines/jakobshavn_xprofile_20090406_points.shp tmp_jakobshavn_flightlines_2008_250m.nc profile_ctr_jakobshavn_flightlines_2008_250m.nc

extract-profile.py ~/data/data_sets/GreenlandFlightlines/jakobshavn_xprofile_20090406_points.shp jakobshavn_flightlines_2008_125m.nc profile_ctr_jakobshavn_flightlines_2008_125m.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_ctr_250m -v thk --labels 'flightline' --figure_title CTR profile_ctr_jakobshavn_flightlines_250m.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_s8_thk_beta_0.8 -v topg --labels '$-$50m,default,$+$50m' --figure_title S8 ../bedmachine_250m_2008_m50m/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_2008_p50m/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc 

~/base/pypismtools/scripts/profile-plot.py -o profile_str_thk_beta_0.8 -v topg --labels '$-$50m,default,$+$50m' --figure_title STR ../bedmachine_250m_2008_m50m/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_2008_p50m/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc 

im-plot.py -p small_font --inner_titles '$-$50m,default,$+$50m' -v topg --singlerow --colormap ~/base/pypismtools/colormaps/wiki-2.0.cpt --colorbar_label -o thk_sensitivity_topg.pdf ../bedmachine_250m_2008_m50m/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_2008_p50m/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc

im-plot.py -p small_font --inner_titles '$-$50m,default,$+$50m' -v divHU --singlerow --colormap RdBu_r  --bounds -500 500 --colorbar_label -o thk_sensitivity_divHU.pdf ../bedmachine_250m_2008_m50m/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_2008_p50m/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc


im-plot.py -p small_font --colorbar_position right --inner_titles '$,$\gamma=0$,$\gamma=0.1$\gamma=0.5$,$\gamma=1$,$\gamma=10$,$\gamma=100$' -v topg --colormap ~/base/pypismtools/colormaps/wiki-2.0.cpt --colorbar_label -o gamma_sensitivity_topg.pdf jakobshavn/jakobshavn_250m_alpha_0.0_gamma_0.0_vscale_0.9_nodhdt_nobmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_0.1_vscale_0.9_nodhdt_nobmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.9_nodhdt_nobmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_nodhdt_nobmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.9_nodhdt_nobmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.9_nodhdt_nobmelt.nc













im-plot.py -p small_font --colorbar_position right --inner_titles '$\gamma=0.5$,$\gamma=1$,$\gamma=5$,$\gamma=10$,$\gamma=0.5$,$\gamma=1$,$\gamma=5$,$\gamma=10$' -v topg --colormap ~/base/pypismtools/colormaps/wiki-2.0.cpt --colorbar_label -o gamma_sensitivity_topg.pdf jakobshavn/jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc   jakobshavn/jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.8_nodhdt_nobmelt.nc 

im-plot.py -p small_font --colorbar_position right --inner_titles '$\gamma=0.5$,$\gamma=1$,$\gamma=5$,$\gamma=10$,$\gamma=0.5$,$\gamma=1$,$\gamma=5$,$\gamma=10$' -v divHU --colormap RdBu_r --bounds -500 500 --colorbar_label -o gamma_sensitivity_divHU.pdf jakobshavn/jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc   jakobshavn/jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.8_nodhdt_nobmelt.nc

im-plot.py -p small_font --inner_titles 'no submarine melt / 1985 data, 900 m/a / 1985 data,no submarine melt / 2008 data, 900 m/a / 2008 data' -v topg --colormap ~/base/pypismtools/colormaps/wiki-2.0.cpt  --colorbar_position right --colorbar_label -o bmelt_topg.pdf ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_dhdt_bmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_dhdt_bmelt.nc

im-plot.py -p small_font --inner_titles 'no submarine melt / 1985 data, 900 m/a / 1985 data,no submarine melt / 2008 data, 900 m/a / 2008 data' -v thk --colormap Accent --bounds 0 2000 --colorbar_position right --colorbar_label -o bmelt_thk.pdf ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_dhdt_bmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_dhdt_bmelt.nc

im-plot.py -p small_font --inner_titles 'no submarine melt / 1985 data, 900 m/a / 1985 data,no submarine melt / 2008 data, 900 m/a / 2008 data' -v divHU --colormap RdBu_r  --bounds -500 500 --colorbar_position right --colorbar_label -o bmelt_divHU.pdf ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_dhdt_bmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_dhdt_bmelt.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_s8_dem_test -v topg --labels 'reference,gimp,cresis' --figure_title S8 profile_s8_jakobshavn_cresis_250m.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_2008_usurf/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_str_dem_test -v topg --labels 'reference,gimp,cresis' --figure_title STR profile_str_jakobshavn_cresis_250m.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_2008_usurf/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc

extract-profile.py ~/data/data_sets/GreenlandFlightlines/jakobshavn_xprofile_20090402-02-08.shp jakobshavn/jakobshavn_250m_alpha_100.1_gamma_0.1_vscale_0.8_dhdt_nobmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_100.1_gamma_0.1_vscale_0.8_dhdt_nobmelt.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_ctr_betas_gamma_1.0  -v topg --labels 'reference,0.7,0.8,0.9,1.0,new' --figure_title CTR  profile_ctr_cresis.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.7_nodhdt_nobmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_nodhdt_nobmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_1.0_nodhdt_nobmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_100.1_gamma_0.1_vscale_0.8_dhdt_nobmelt.nc


~/base/pypismtools/scripts/profile-plot.py -o profile_ctr_betas_gamma_1.0  -v thk --labels 'reference,0.7,0.8,0.9,1.0,new' --figure_title CTR  profile_ctr_cresis.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.7_nodhdt_nobmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.9_nodhdt_nobmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_1.0_nodhdt_nobmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_100.1_gamma_0.1_vscale_0.8_dhdt_nobmelt.nc

im-plot.py -p twocol --colorbar_label -v topg --inner_titles '0,CReSIS' --singlerow --colormap ~/base/pypismtools/colormaps/wiki-2.0.cpt -o topg_cresis_0_topg_250m.pdf jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_2008/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc

im-plot.py -p twocol --colorbar_label -v thk  --bounds -200 200 --colormap RdBu_r --singlerow -o rho_diff_thk_250m.pdf --obs_file withrho_gamma_1.0.nc norho_gamma_1.0.nc

im-plot.py -p twocol --colorbar_label -v thk  --bounds -100 100 --colormap RdBu_r --singlerow -o diff.pdf --obs_file jakobshavn/jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_0.5_vscale_0.8_dhdt_bmelt.nc 
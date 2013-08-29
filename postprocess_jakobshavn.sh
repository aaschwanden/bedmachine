#!/bin/bash

set -x

R=l
GS=250
fill_value=-2e9

alpha=0.0

for gamma in 0.0 0.5 1.0 2.0 5.0 10.0 20.0 50.0 100.0 200.0 500.0 1000.0 2000.0 5000.0 10000.0
do
    for project in "jakobshavn"
    do
        for exp in  "nodhdt_nobmelt" "dhdt_nobmelt"
        do
	    for scale in 0.2 0.4 0.6 0.8 1.0
	    do
		nc2cdo.py ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc
		extract-profile.py ~/data/data_sets/Jakobshavn1985/jakobshavn_xprofile_str_point.shp ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc ${project}/profile_str_${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc
		extract-profile.py ~/data/data_sets/Jakobshavn1985/jakobshavn_xprofile_s8_point.shp ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc ${project}/profile_s8_${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc
		extract-profile.py ~/data/data_sets/GreenlandFlightlines/jakobshavn_xprofile_20090406_points.shp ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc ${project}/profile_ctr_${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc
	    done
	done
    done
done

exit

for gamma in 0.0 0.5 1.0 2.0 5.0 10.0 20.0 50.0 100.0 200.0 500.0 1000.0 2000.0 5000.0 10000.0
do
    for project in "jakobshavn"
    do
        for exp in  "nodhdt_nobmelt" "dhdt_nobmelt"
        do
	    for scale in 0.2 0.4 0.6 0.8 1.0
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

for project in "jakobshavn"
do
    for dataset in "cresis" "umt" "flightlines"
    do
        # thk of input data sets
        #im-plot.py -p medium --colorbar_label --bounds 0 2500 -v thk --singlerow --colormap /Users/andy/base/PyPISMTools/colormaps/sst.cpt -o ${project}/${project}_${dataset}_${GS}m_thk.pdf ${project}_${dataset}_${GS}m.nc
        # topg of input data sets
        #im-plot.py -p medium --colorbar_label -v topg --singlerow --colormap /Users/andy/base/pypismtools/colormaps/wiki-2.0.cpt -o ${project}/${project}_${dataset}_${GS}m_topg.pdf ${project}_${dataset}_${GS}m.nc
        extract-profile.py ~/data/data_sets/Jakobshavn1985/jakobshavn_xprofile_str_point.shp ${project}_${dataset}_${GS}m.nc profile_str_${project}_${dataset}_${GS}m.nc
        extract-profile.py ~/data/data_sets/Jakobshavn1985/jakobshavn_xprofile_s8_point.shp ${project}_${dataset}_${GS}m.nc profile_s8_${project}_${dataset}_${GS}m.nc
        extract-profile.py ~/data/data_sets/GreenlandFlightlines/jakobshavn_xprofile_20090406_points.shp ${project}_${dataset}_${GS}m.nc profile_ctr_${project}_${dataset}_${GS}m.nc
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


~/base/pypismtools/scripts/profile-plot.py -o profile_s8_gammas_beta_0.8 --split 2 -v topg --labels 'reference,1,2,5,10,100,1000' --figure_title S8 profile_s8_jakobshavn_cresis_250m.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_2.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1000.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/profile_s8_jakobshavn_cresis_250m.nc ../bedmachine_250m_1985/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_2.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1000.0_vscale_0.8_nodhdt_nobmelt.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_str_gammas_beta_0.8 --split 2 -v topg --labels 'reference,1,2,5,10,100,1000' --figure_title STR profile_str_jakobshavn_cresis_250m.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_2.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1000.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/profile_str_jakobshavn_cresis_250m.nc ../bedmachine_250m_1985/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_2.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1000.0_vscale_0.8_nodhdt_nobmelt.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_ctr_gammas_beta_0.8 -v thk --labels 'flightline,reference,1,5,10,100,1000' --figure_title CTR profile_ctr_jakobshavn_flightlines_2008_250m.nc profile_ctr_jakobshavn_cresis_250m.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_1000.0_vscale_0.8_nodhdt_nobmelt.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_ctr_gammas_beta_0.8 --split 2 -v topg --labels 'reference,1,2,5,10,100,1000' --figure_title CTR profile_ctr_jakobshavn_cresis_250m.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_2.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_1000.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/profile_ctr_jakobshavn_cresis_250m.nc ../bedmachine_250m_1985/jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_2.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_ctr_jakobshavn_250m_alpha_0.0_gamma_1000.0_vscale_0.8_nodhdt_nobmelt.nc


~/base/pypismtools/scripts/profile-plot.py -o profile_s8_betas_gamma_5.0 --split 2 -v topg --labels 'reference,0.6,0.8,1.0' --figure_title S8 profile_s8_jakobshavn_cresis_250m.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.6_nodhdt_nobmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_1.0_nodhdt_nobmelt.nc profile_s8_jakobshavn_cresis_250m.nc ../bedmachine_250m_1985/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.6_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_1.0_nodhdt_nobmelt.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_str_betas_gamma_5.0 --split 2 -v topg --labels 'reference,0.6,0.8,1.0' --figure_title STR profile_str_jakobshavn_cresis_250m.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.6_nodhdt_nobmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_1.0_nodhdt_nobmelt.nc profile_str_jakobshavn_cresis_250m.nc ../bedmachine_250m_1985/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.6_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_1.0_nodhdt_nobmelt.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_s8_betas_gamma_1.0 --split 2 -v topg --labels 'reference,0.6,0.8,1.0' --figure_title S8 profile_s8_jakobshavn_cresis_250m.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.6_nodhdt_nobmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_1.0_nodhdt_nobmelt.nc profile_s8_jakobshavn_cresis_250m.nc ../bedmachine_250m_1985/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.6_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_1.0_nodhdt_nobmelt.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_str_betas_gamma_1.0 --split 2 -v topg --labels 'reference,0.6,0.8,1.0' --figure_title STR profile_str_jakobshavn_cresis_250m.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.6_nodhdt_nobmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_1.0_nodhdt_nobmelt.nc profile_str_jakobshavn_cresis_250m.nc ../bedmachine_250m_1985/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.6_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_1.0_nodhdt_nobmelt.nc


im-plot.py -p twocol --colorbar_label -v magnitude --inner_titles '1985,2008' --singlerow --colormap /Users/andy/base/pypismtools/colormaps/Full_saturation_spectrum_CCW_orange.cpt -o speed_1985_2008.pdf jakobshavn_surf_vels_1985_250m.nc jakobshavn_surf_vels_2008_2009_250m.nc

im-plot.py --bounds -1000 2100 -p twocol --colorbar_label -v usurf --colormap ~/base/pypismtools/colormaps/wiki-2.0.cpt --inner_titles '1985,2008' --singlerow -o usurf_1985_2008.pdf jakobshavn_usurf_1985_250m.nc jakobshavn_usurf_2008_250m.nc

im-plot.py -p twocol --colorbar_label -v topg --inner_titles 'Bamber 2001,CReSIS' --singlerow --colormap ~/base/pypismtools/colormaps/wiki-2.0.cpt -o topg_bamber_cresis_topg_250m.pdf jakobshavn_searise_v1.1_250m.nc jakobshavn_cresis_250m.nc

gdalwarp -overwrite -r bilinear -te -220000 -2300000 -120000 -2240000 -tr 100 100 -t_srs EPSG:3413 /Volumes/data/data_sets/GreenlandLandsat/LC80090112013101LGN01_pansharpended_natural_z.tif kansas/jakobshavn_landsat8.tif

basemap-plot.py --geotiff_file jakobshavn_gimpdem_90m_hillshade.tif  -p twocol --colorbar_label -v topg --inner_title --singlerow --colormap ~/base/pypismtools/colormaps/wiki-2.0.cpt -o topg_bamber_cresis_topg_250m.pdf Greenland_5km_v1.1.nc jakobshavn_cresis_250m.nc

~/base/pypismtools/scripts/basemap-plot.py --geotiff_file jakobshavn_gimpdem_90m_hillshade.tif --alpha 0.5 -p twocol --colorbar_label -v topg --inner_titles 'Bamber (2001),Bamber (2013),CReSIS (2012)' --singlecol --colormap ~/base/pypismtools/colormaps/wiki-2.0.cpt -o bamber_cresis_topg_250m.pdf Greenland_5km_v1.1.nc ~/data/data_sets/greenland_beds_v2/Greenland_1km_v2.nc jakobshavn_cresis_250m.nc

project=jakobshavn
GS=250
alpha=0.0
exp=nodhdt_nobmelt
for scale in 0.2 0.4 0.6 0.8 1.0;
do
    im-plot.py -p twocol --colorbar_position right --inner_titles '$\gamma$=1,$\gamma$=2,$\gamma$=5,$\gamma$=10,$\gamma$=100,$\gamma$=1000' --colorbar_label -v divHU --colormap RdBu_r  --bounds -500 500 -o ${project}_${GS}m_alpha_${alpha}_gammas_vscale_${scale}_${exp}_divHU.pdf ${project}/${project}_${GS}m_alpha_${alpha}_gamma_1.0_vscale_${scale}_${exp}.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_2.0_vscale_${scale}_${exp}.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_5.0_vscale_${scale}_${exp}.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_10.0_vscale_${scale}_${exp}.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_100.0_vscale_${scale}_${exp}.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_1000.0_vscale_${scale}_${exp}.nc
done


extract-profile.py ~/data/data_sets/GreenlandFlightlines/jakobshavn_xprofile_20090406_points.shp tmp_jakobshavn_flightlines_2008_250m.nc profile_ctr_jakobshavn_flightlines_2008_250m.nc

extract-profile.py ~/data/data_sets/GreenlandFlightlines/jakobshavn_xprofile_20090406_points.shp jakobshavn_flightlines_2008_125m.nc profile_ctr_jakobshavn_flightlines_2008_125m.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_ctr_250m -v thk --labels 'flightline' --figure_title CTR profile_ctr_jakobshavn_flightlines_250m.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_s8_thk_beta_0.8 -v topg --labels '$-$50m,default,$+$50m' --figure_title S8 ../bedmachine_250m_2008_m50m/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_2008_p50m/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc 

~/base/pypismtools/scripts/profile-plot.py -o profile_str_thk_beta_0.8 -v topg --labels '$-$50m,default,$+$50m' --figure_title STR ../bedmachine_250m_2008_m50m/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_2008_p50m/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc 

im-plot.py -p small_font --inner_titles '$-$50m,default,$+$50m' -v topg --singlerow --colormap ~/base/pypismtools/colormaps/wiki-2.0.cpt --colorbar_label -o thk_sensitivity_topg.pdf ../bedmachine_250m_2008_m50m/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_2008_p50m/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc

im-plot.py -p small_font --inner_titles '$-$50m,default,$+$50m' -v divHU --singlerow --colormap RdBu_r  --bounds -500 500 --colorbar_label -o thk_sensitivity_divHU.pdf ../bedmachine_250m_2008_m50m/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_2008_p50m/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc

im-plot.py -p small_font --colorbar_position right --inner_titles '$\gamma=1$,$\gamma=2$,$\gamma=5$,$\gamma=10$,$\gamma=1$,$\gamma=2$,$\gamma=5$,$\gamma=10$' -v topg --colormap ~/base/pypismtools/colormaps/wiki-2.0.cpt --colorbar_label -o gamma_sensitivity_topg.pdf jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc  jakobshavn/jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.8_nodhdt_nobmelt.nc

im-plot.py -p small_font --colorbar_position right --inner_titles '$\gamma=1$,$\gamma=2$,$\gamma=5$,$\gamma=10$,$\gamma=1$,$\gamma=2$,$\gamma=5$,$\gamma=10$' -v divHU --colormap RdBu_r --bounds -500 500 --colorbar_label -o gamma_sensitivity_divHU.pdf jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc  jakobshavn/jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.8_nodhdt_nobmelt.nc

im-plot.py -p small_font --inner_titles 'no submarine melt / 1985 data, 600 m/a / 1985 data,no submarine melt / 2008 data, 600 m/a / 2008 data' -v topg --colormap ~/base/pypismtools/colormaps/wiki-2.0.cpt  --colorbar_position right --colorbar_label -o bmelt_topg.pdf ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_2.0_gamma_1.0_vscale_0.8_nodhdt_bmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/jakobshavn_250m_alpha_2.0_gamma_1.0_vscale_0.8_nodhdt_bmelt.nc

im-plot.py -p small_font --inner_titles 'no submarine melt / 1985 data, 600 m/a / 1985 data,no submarine melt / 2008 data, 600 m/a / 2008 data' -v thk --colormap Accent --bounds 0 2000 --colorbar_position right --colorbar_label -o bmelt_thk.pdf ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/jakobshavn_250m_alpha_2.0_gamma_1.0_vscale_0.8_nodhdt_bmelt.nc jakobshavn/jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc jakobshavn/jakobshavn_250m_alpha_2.0_gamma_1.0_vscale_0.8_nodhdt_bmelt.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_s8_dem_test -v topg --labels 'reference,gimp,cresis' --figure_title S8 profile_s8_jakobshavn_cresis_250m.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_2008_usurf/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_str_dem_test -v topg --labels 'reference,gimp,cresis' --figure_title STR profile_str_jakobshavn_cresis_250m.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc ../bedmachine_250m_2008_usurf/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.8_nodhdt_nobmelt.nc
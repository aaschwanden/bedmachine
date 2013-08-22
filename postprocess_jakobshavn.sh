#!/bin/bash

set -x

R=l
GS=250
fill_value=-2e9

alpha=0.0

for gamma in 0.0 1.0 2.0 5.0 10.0 20.0 50.0 100.0 200.0 500.0 1000.0 2000.0 5000.0 10000.0 20000.0 50000.0
do
    for project in "jakobshavn"
    do
        for exp in  "nodhdt_nobmelt" "dhdt_nobmelt"
        do
	    for scale in 0.25 0.5 0.75 1.0 1.25
	    do
		nc2cdo.py ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc
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
		extract-profile.py jakobshavn_xprofile_str_point.shp ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc ${project}/profile_str_${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc
		extract-profile.py jakobshavn_xprofile_s8_point.shp ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc ${project}/profile_s8_${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc
	    done
	done
    done
done

for project in "jakobshavn"
do
    for dataset in "cresis" "umt"
    do
        # thk of input data sets
        im-plot.py -p medium --colorbar_label --bounds 0 2500 -v thk --singlerow --colormap /Users/andy/base/PyPISMTools/colormaps/sst.cpt -o ${project}/${project}_${dataset}_${GS}m_thk.pdf ${project}_${dataset}_${GS}m.nc
        # topg of input data sets
        im-plot.py -p medium --colorbar_label -v topg --singlerow --colormap /Users/andy/base/pypismtools/colormaps/wiki-2.0.cpt -o ${project}/${project}_${dataset}_${GS}m_topg.pdf ${project}_${dataset}_${GS}m.nc
        extract-profile.py jakobshavn_xprofile_str_point.shp ${project}_${dataset}_${GS}m.nc profile_str_${project}_${dataset}_${GS}m.nc
        extract-profile.py jakobshavn_xprofile_s8_point.shp ${project}_${dataset}_${GS}m.nc profile_s8_${project}_${dataset}_${GS}m.nc
    done
done

scale=0.75
alpha=0.0
exp=nodhdt_nobmelt
project=jakobshavn
gamma=5.0
GS=250
im-plot.py -p twocol --inner_title --colorbar_label -v divHU_cresis --singlerow --colormap RdBu_r  --bounds -500 500 -o ${project}_${GS}m_cresis_vscale_${scale}_${exp}_divHU.pdf ../bedmachine_250m_1985/${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc 


 im-plot.py -p medium --colorbar_label -v divHU_cresis --singlerow --colormap RdBu_r  --bounds -500 500 -o ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}_divHU_cresis.pdf ${project}/${project}_${GS}m_alpha_${alpha}_gamma_${gamma}_vscale_${scale}_${exp}.nc

extract-profile.py jakobshavn_xprofile_str_point.shp ${project}_${dataset}_${GS}m.nc profile_str_${project}_${dataset}_${GS}m.nc


~/base/pypismtools/scripts/profile-plot.py -o profile_s8_gamma --split 2 -v topg --labels 'reference,1,2,5,10,100,1000' --figure_title S8 profile_s8_jakobshavn_cresis_250m.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.75_nodhdt_nobmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_2.0_vscale_0.75_nodhdt_nobmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.75_nodhdt_nobmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.75_nodhdt_nobmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.75_nodhdt_nobmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1000.0_vscale_0.75_nodhdt_nobmelt.nc ../bedmachine_250m_1985/profile_s8_jakobshavn_cresis_250m.nc ../bedmachine_250m_1985/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.75_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_2.0_vscale_0.75_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.75_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.75_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.75_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_1000.0_vscale_0.75_nodhdt_nobmelt.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_str_gamma --split 2 -v topg --labels 'reference,1,2,5,10,100,1000' --figure_title STR profile_str_jakobshavn_cresis_250m.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.75_nodhdt_nobmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_2.0_vscale_0.75_nodhdt_nobmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.75_nodhdt_nobmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.75_nodhdt_nobmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.75_nodhdt_nobmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1000.0_vscale_0.75_nodhdt_nobmelt.nc ../bedmachine_250m_1985/profile_str_jakobshavn_cresis_250m.nc ../bedmachine_250m_1985/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1.0_vscale_0.75_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_2.0_vscale_0.75_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.75_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_10.0_vscale_0.75_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_100.0_vscale_0.75_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_1000.0_vscale_0.75_nodhdt_nobmelt.nc


~/base/pypismtools/scripts/profile-plot.py -o profile_s8_alpha --split 2 -v topg --labels 'reference,0.25,0.5,0.75,1.0,1.25' --figure_title S8 profile_s8_jakobshavn_cresis_250m.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.25_nodhdt_nobmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.5_nodhdt_nobmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.75_nodhdt_nobmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_1.0_nodhdt_nobmelt.nc jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_1.25_nodhdt_nobmelt.nc profile_s8_jakobshavn_cresis_250m.nc ../bedmachine_250m_1985/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.25_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.5_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.75_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_1.0_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_s8_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_1.25_nodhdt_nobmelt.nc

~/base/pypismtools/scripts/profile-plot.py -o profile_str_alpha --split 2 -v topg --labels 'reference,0.25,0.5,0.75,1.0,1.25' --figure_title STR profile_str_jakobshavn_cresis_250m.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.25_nodhdt_nobmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.5_nodhdt_nobmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.75_nodhdt_nobmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_1.0_nodhdt_nobmelt.nc jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_1.25_nodhdt_nobmelt.nc profile_str_jakobshavn_cresis_250m.nc ../bedmachine_250m_1985/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.25_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.5_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_0.75_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_1.0_nodhdt_nobmelt.nc ../bedmachine_250m_1985/jakobshavn/profile_str_jakobshavn_250m_alpha_0.0_gamma_5.0_vscale_1.25_nodhdt_nobmelt.nc


im-plot.py -p twocol --colorbar_label -v magnitude --inner_title --singlerow --colormap /Users/andy/base/pypismtools/colormaps/Full_saturation_spectrum_CCW_orange.cpt -o speed_1985_2008.pdf jakobshavn_surf_vels_1985_250m.nc jakobshavn_surf_vels_2008_2009_250m.nc

im-plot.py -p twocol --colorbar_label -v usurf --inner_title --singlerow -o usurf_1985_2008.pdf jakobshavn_usurf_1985_250m.nc jakobshavn_usurf_2008_250m.nc

im-plot.py -p twocol --colorbar_label -v topg --inner_title --singlerow --colormap ~/base/pypismtools/colormaps/wiki-2.0.cpt -o topg_bamber_cresis_topg_250m.pdf jakobshavn_searise_v1.1_250m.nc jakobshavn_cresis_250m.nc

gdalwarp -overwrite -r bilinear -te -220000 -2300000 -120000 -2240000 -tr 100 100 -t_srs EPSG:3413 /Volumes/data/data_sets/GreenlandLandsat/LC80090112013101LGN01_pansharpended_natural_z.tif kansas/jakobshavn_landsat8.tif

basemap-plot.py --geotiff_file jakobshavn_gimpdem_90m_hillshade.tif  -p twocol --colorbar_label -v topg --inner_title --singlerow --colormap ~/base/pypismtools/colormaps/wiki-2.0.cpt -o topg_bamber_cresis_topg_250m.pdf Greenland_5km_v1.1.nc jakobshavn_cresis_250m.nc

~/base/pypismtools/scripts/basemap-plot.py --geotiff_file jakobshavn_gimpdem_90m_hillshade.tif  -p twocol --colorbar_label -v topg --inner_title --singlecol --colormap ~/base/pypismtools/colormaps/wiki-2.0.cpt -o bamber_cresis_topg_250m.pdf Greenland_5km_v1.1.nc jakobshavn_Greenland_1km_v2.nc jakobshavn_cresis_250m.nc
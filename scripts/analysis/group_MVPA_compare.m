%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to compare results from different MVPA
% searchlight analyses in the same coordinate space
% Tom Possidente - July 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
addpath(genpath('/projectnb/somerslab/tom/ArthurfMRI-main/'));
ccc;

%% Set up key variables and paths

basedir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/MVPA_results/';
reg_file = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/fsaverage/mri.2mm/reg.2mm.dat'; % needed to transform from vol to surf

% contrast1_masked = 'tstats_mean_masked_VWM_VSMC+AWM_ASMC';
% contrast2_masked = 'tstats_mean_masked_VSMC_F+ASMC_F';
% contrast1_unmasked = 'tstats_mean_VWM_VSMC+AWM_ASMC';
% contrast2_unmasked = 'tstats_mean_VSMC_F+ASMC_F';
% pval_thresh = 0.9999; 
contrast1_masked = 'Vis_F';
contrast2_masked = 'Aud_F';
contrast1_unmasked = 'Vis_F';
contrast2_unmasked = 'Aud_F';
pval_thresh = 0.05; 

hemispheres = {'lh', 'rh'};

%% Load maps, take union, take difference
for hh = 1:length(hemispheres)
    hemi = hemispheres{hh};

    %cont1_pmap = MRIread([basedir contrast1_masked '_blocktrials_MVPA_FWE05_TFCE_surf_' hemi '.nii']);
    cont1_pmap = MRIread([basedir contrast1_masked '_blocktrials_MVPA_FWE_tfce_corrp_tstat1_surf_' hemi '.nii']);
    cont1_pmap.vol = cont1_pmap.vol >= 1-pval_thresh; % inv pmap
    %cont2_pmap = MRIread([basedir contrast2_masked '_blocktrials_MVPA_FWE05_TFCE_surf_' hemi '.nii']);
    cont2_pmap = MRIread([basedir contrast2_masked '_blocktrials_MVPA_FWE_tfce_corrp_tstat1_surf_' hemi '.nii']);
    cont2_pmap.vol = cont2_pmap.vol >= 1-pval_thresh; % inv pmap
    
    pmap_union = cont1_pmap.vol | cont2_pmap.vol; % OR operation on pmaps

    % Load tstat maps
    % cont1_tmap = MRIread([basedir contrast1_unmasked '_blocktrials_MVPA_FWE05_TFCE_surf_' hemi '.nii']);
    % cont2_tmap = MRIread([basedir contrast2_unmasked '_blocktrials_MVPA_FWE05_TFCE_surf_' hemi '.nii']);
    cont1_tmap = MRIread([basedir contrast1_unmasked '_blocktrials_MVPA_FWE_tstat1_surf_' hemi '.nii']);
    cont2_tmap = MRIread([basedir contrast2_unmasked '_blocktrials_MVPA_FWE_tstat1_surf_' hemi '.nii']);
    tmaps_neg = cont1_tmap.vol < 0 | cont2_tmap.vol < 0;

    tmap_diff = cont1_tmap.vol - cont2_tmap.vol; % take diff of tstat maps
    tmap_diff(~pmap_union | tmaps_neg) = 0; % make all tstats 0 for areas not significant in either maps and areas with negative tstats in either map
    
    tmap_diff_signmasked = MRIread('/share/pkg.8/freesurfer/7.4.1_CentOS-8/install/freesurfer/subjects/fsaverage/mri.2mm/mni305.cor.mgz'); % load template nii
    tmap_diff_signmasked.vol = tmap_diff;
    % MRIwrite(tmap_diff_signmasked, [basedir 'tstats_diff_VWM_VSMC+AWM_ASMC-VSMC_F+ASMC_F_blocktrials_MVPA_FWE05_TFCE_' hemi '.nii']);
    MRIwrite(tmap_diff_signmasked, [basedir 'tstats_diff_' contrast1_masked '-' contrast2_masked '_blocktrials_MVPA_FWE05_TFCE_' hemi '.nii']);

end
    


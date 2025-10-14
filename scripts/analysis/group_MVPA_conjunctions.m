%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to make conjunction maps from the MVPA
% searchlight results
% Tom Possidente - August 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
addpath(genpath('/projectnb/somerslab/tom/ArthurfMRI-main/'));
ccc;

%% Set up key variables and paths

basedir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/MVPA_results/';
reg_file = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/fsaverage/mri.2mm/reg.2mm.dat'; % needed to transform from vol to surf

contrast1 = 'crosstask_WM-SMC_SMC-F';
contrast2 = 'crosstask_SMC-F_WM-SMC';
pval_thresh = 0.05; 

hemispheres = {'lh', 'rh'};

%% Load MVPA pmaps and conjoin 
for hh = 1:length(hemispheres)
    hemi = hemispheres{hh};

    cont1_pmap = MRIread([basedir contrast1 '_blocktrials_MVPA_FWE_tfce_corrp_tstat1_surf_' hemi '.nii']);
    cont1_pmap.vol = cont1_pmap.vol >= 1-pval_thresh; % inv pmap
    cont2_pmap = MRIread([basedir contrast2 '_blocktrials_MVPA_FWE_tfce_corrp_tstat1_surf_' hemi '.nii']);
    cont2_pmap.vol = cont2_pmap.vol >= 1-pval_thresh; % inv pmap
    
    pmap_combined = cont1_pmap.vol & cont2_pmap.vol; % AND operation on pmaps

    % Load tstat maps
    cont1_tmap = MRIread([basedir contrast1 '_blocktrials_MVPA_FWE_tstat1_surf_' hemi '.nii']);
    cont2_tmap = MRIread([basedir contrast2 '_blocktrials_MVPA_FWE_tstat1_surf_' hemi '.nii']);
    tmap_mean = mean([cont1_tmap.vol; cont2_tmap.vol]); % take mean of tstat maps

    % Save out mean of tstat maps unmasked
    tmap_mean_unmasked = MRIread('/share/pkg.8/freesurfer/7.4.1_CentOS-8/install/freesurfer/subjects/fsaverage/mri.2mm/mni305.cor.mgz'); % load template nii
    tmap_mean_unmasked.vol = tmap_mean;
    MRIwrite(tmap_mean_unmasked, [basedir 'tstats_mean_' contrast1 '+' contrast2 '_blocktrials_MVPA_FWE05_TFCE_surf_' hemi '.nii']);

    % make all tstats 0 for areas not significant in both maps and save out
    tmap_mean(~pmap_combined) = 0; 
    tmap_mean_sigmasked = MRIread('/share/pkg.8/freesurfer/7.4.1_CentOS-8/install/freesurfer/subjects/fsaverage/mri.2mm/mni305.cor.mgz'); % load template nii
    tmap_mean_sigmasked.vol = tmap_mean;
    MRIwrite(tmap_mean_sigmasked, [basedir 'tstats_mean_masked_' contrast1 '+' contrast2 '_blocktrials_MVPA_FWE05_TFCE_surf_' hemi '.nii']);
end
    


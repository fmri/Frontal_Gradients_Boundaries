%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to perform searchlight MVPA analysis on the
% trifloc localizer data to predict presence of each condition 
% Tom Possidente - July 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
addpath(genpath('/projectnb/somerslab/tom/ArthurfMRI-main/'));
ccc;

%% Set up key variables and paths
% Load subj codes
experiment_name = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode(~ismember(subjDf_cut.subjCode, {'AH', 'SL', 'RR', 'PP', 'MM'})); % N=19 (all have fixation condition)
%subjCodes = subjDf_cut.subjCode(~ismember(subjDf_cut.subjCode, {'AH', 'SL', 'RR'})); % N=21 (MM and PP don't have fixation)
N_subjs = length(subjCodes);

basedir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';
results_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/MVPA_results/';

% We are trying to predict WM vs SMC trials, so we are looking at all WM and all SMC conditions
conds_perblock_str = 'blocktrials'; % blocktrials or 16trials
save_str = 'Aud_F';
conds_perblock = 1; %16;
conditions = [4,5,6,10]; %[1:96]; % all visual WM, auditory WM, visual SMC, and auditory SMC
condition_labels = [1,1,1,0]; %[ones(conds_perblock*2,1); zeros(conds_perblock,1); ones(conds_perblock*2,1); zeros(conds_perblock,1)];

% localizer condition order (sections_per_block=2):
% 1-4 = vA (block 1, section 1; block 1, section 2; block 2, section 1; block 2, section 2)
% 5-6 = vP (block 1, section 1; block 1, section 2)
% 7-10 = aA (block 1, section 1; block 1, section 2; block 2, section 1; block 2, section 2)
% 11-12 = aP (block 1, section 1; block 1, section 2)
% 13-16 = tA (block 1, section 1; block 1, section 2; block 2, section 1; block 2, section 2)
% 17-18 = tP (block 1, section 1; block 1, section 2)
% 19-20 = f (block 1, section 1; block 1, section 2)
% localizer condition order (sections_per_block=16):
% 1-32 = vA 
% 33-48 = vP 
% 49-80 = aA 
% 81-96 = aP 
% 97-128 = tA 
% 129-144 = tP 
% 145-160 = f 

if any(ismember(conditions, 10)) && any(ismember(subjCodes, {'MM', 'PP'}))
    error('MM and PP are included in subjCodes and fixation is included in this analysis, but MM and PP do not have fixation conditions. Remove MM and PP from subjCodes');
end
if ismember('PP', subjCodes)
    trials_mod = 2; % this is the only subj with 6 runs, so use mod of 2 (# of extra runs)
else
    trials_mod = 0;
end

% Load cerebral cortex mask
cortex_mask = MRIread('/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/fsaverage/mri.2mm/grey_cortex_mask_dilate2erode1.2mm.mgz');
cortex_mask = cortex_mask.vol>0;
cortex_mask_reshaped = cortex_mask(:);

% Set up storage variables
N_trials_perrun = length(conditions); 
N_runs = 4; % all subjs have 4 runs (except PP has 6, this will be accounted for)
N_trials_persubj = N_trials_perrun * N_runs; 
x = nan((N_subjs * N_trials_persubj)+(trials_mod*N_trials_perrun), sum(cortex_mask_reshaped,'all'));
y = nan((N_subjs * N_trials_persubj)+(trials_mod*N_trials_perrun), 1);
subjlist = nan((N_subjs * N_trials_persubj)+(trials_mod*N_trials_perrun), 1);
runlist = nan((N_subjs * N_trials_persubj)+(trials_mod*N_trials_perrun), 1); % I don't think we will actually need this, but just in case

%% For each subj/run, load GLM coefficient data and z-score
count = 0;
for ss = 1:N_subjs
    subjCode = subjCodes{ss};
    N_runs_mod = 0;

    if ismember(subjCode, {'MM', 'PP'})
        suffix = ['_' num2str(conds_perblock) 'section_nofix'];
        if strcmp(subjCode,'PP')
            N_runs_mod = 2;
        end
    else
        suffix = ['_' num2str(conds_perblock) 'section'];
    end

    for rr = 1:N_runs+N_runs_mod
        count = count + 1;

        if strcmp(subjCode, 'LA') && rr == 2 && strcmp(conds_perblock_str, '16trials')
            continue;
        end

        % Load GLM coefficients
        coefs = MRIread([basedir subjCode '/localizer/localizer_block_regression_0sm_run' num2str(rr) suffix '/beta.nii.gz']);

        % Load param file to get coef order
        %conds = readmatrix([basedir subjCode '/localizer/00' num2str(rr) '/localizer_condition_timing_perblock_1section.para'], 'FileType','text');

        % Get coefs for relevant conditions
        coefs = coefs.vol(:,:,:,conditions);
        coefs_reshaped = reshape(coefs, size(coefs,1)*size(coefs,2)*size(coefs,3), length(conditions));
        coefs_inmask = coefs_reshaped(cortex_mask_reshaped==1,:);

        % Z-score coefs
        coefs_zscored = zscore(coefs_inmask'); % This is z-scored each column (voxel) separately

        % Store data
        index_start = ((count-1)*N_trials_perrun + 1);
        indices = index_start:(index_start+N_trials_perrun-1);
        x(indices,:) = coefs_zscored;
        y(indices) = condition_labels; % always the same because coefs are stored in the same condition order (1 = WM, 0 = SMC)
        subjlist(indices) = repmat(ss, N_trials_perrun,1);
        runlist(indices) = repmat(rr, N_trials_perrun,1);

    end
    disp(['Finished subj ' subjCode])
end

%save('trial_coefs_WM_SMC_pertrial.mat', 'x', 'y', 'subjlist', 'runlist');

%load('/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/MVPA_results/trial_coefs_SMC_F.mat', 'x', 'y', 'subjlist', 'runlist');

%% Run searchlight analysis
CVfold = subjlist; % information for leave-one-subject-out
dist = 2; % radius of searchlight, in voxels. radius of 2 gives a small sphere that is 33 voxels
metric = 'Pearson'; % will use correlation between choice and prediction to gauge performance
method = 'PLS1'; 

% Remove nans (this should only happen if a full run is bad)
if any(isnan(x),'all')
    nanmask = any(isnan(x),2);
    check_nans = sum(nanmask);
    assert(mod(check_nans, N_trials_perrun)==0, 'sum of nans are not a multiple of trials per run');
    x = x(~nanmask,:);
    y = y(~nanmask);
    CVfold = CVfold(~nanmask);
end

tic;
measure = Afmri.searchlight(x,y,CVfold,cortex_mask,dist,metric,method); % nested leave-one-subject-out cross-validation to find searchlights that can predict Y
toc


% Mask out voxels with any nans
nan_measure_mask = ~any(isnan(measure),1);
measure = measure(:,nan_measure_mask);
cortex_mask_reshaped(cortex_mask_reshaped) = nan_measure_mask;
cortex_mask = reshape(cortex_mask_reshaped, size(cortex_mask,1), size(cortex_mask,2), size(cortex_mask,3));

% measure is N_subj brain images of pearson correlation coefficients, each of which has been predicted out-of-sample by a predictor trained on n-1 subjects
% since correlations are null at zero, we can use sign-flipping permutation test to identify brain regions that can significantly decode

%% Save out files for permutation testing
template = MRIread('/share/pkg.8/freesurfer/7.4.1_CentOS-8/install/freesurfer/subjects/fsaverage/mri.2mm/mni305.cor.mgz');

measure_data = template;
measure_data.vol = zeros(size(measure_data.vol, 1), size(measure_data.vol, 2), size(measure_data.vol, 3), size(measure,1));
measure_data.nframes = size(measure,1);
cortex_mask_rep = repmat(cortex_mask, 1, 1, 1, size(measure,1));
measure_data.vol(cortex_mask_rep) = measure';
MRIwrite(measure_data, [results_dir 'corr_measure_' save_str '_' conds_perblock_str '.nii']);

mask_data = template;
mask_data.vol(cortex_mask) = 1;
mask_data.vol(~cortex_mask) = 0;
MRIwrite(mask_data, [results_dir 'mask_' save_str '_' conds_perblock_str '.nii']);

%% Run TFCE permutation testing
output_file = [results_dir save_str '_' conds_perblock_str '_MVPA_FWE'];
input_data = [results_dir 'corr_measure_' save_str '_' conds_perblock_str '.nii'];
mask = [results_dir 'mask_' save_str '_' conds_perblock_str '.nii'];
disp(['randomise -i ' input_data ' -o ' output_file ' -m ' mask ' -1 -T -n 1000'])
tic;
unix(['randomise -i ' input_data ' -o ' output_file ' -m ' mask ' -1 -T -n 1000']) % using randomise TFCE from FSL with -1 for 1 sample t-test, -T for TFCE, and -n 5000 for 5000 iterations permutation testing
toc

% Convert volume to surface with vol2sur
reg_file = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/fsaverage/mri.2mm/reg.2mm.dat';
unix(['mri_vol2surf --src ' output_file '_tfce_corrp_tstat1.nii.gz --out ' output_file '_tfce_corrp_tstat1_surf_lh.nii --reg ' reg_file ' --trgsubject fsaverage --hemi lh' ])
unix(['mri_vol2surf --src ' output_file '_tfce_corrp_tstat1.nii.gz --out ' output_file '_tfce_corrp_tstat1_surf_rh.nii --reg ' reg_file ' --trgsubject fsaverage --hemi rh' ])
unix(['mri_vol2surf --src ' output_file '_tstat1.nii.gz --out ' output_file '_tstat1_surf_lh.nii --reg ' reg_file ' --trgsubject fsaverage --hemi lh' ])
unix(['mri_vol2surf --src ' output_file '_tstat1.nii.gz --out ' output_file '_tstat1_surf_rh.nii --reg ' reg_file ' --trgsubject fsaverage --hemi rh' ])







%% Old permutation testing
% initp = .001; % a more reasonable threshold, since MVPA is more powerful than GLM
% covar = []; % no covariates to declare. a simple t-test
% stattype = 'size'; % cluster-size testing
% correctiontype = 'FDR'; % False discovery rate correction
% corrtype = []; % doesn't apply here...
% filename = 'group_WM_SMC_decode_searchlight_16trials'; 
% 
% tic
% Afmri.permTest(measure,cortex_mask,initp,covar,stattype,correctiontype,corrtype,filename);
% toc;
% 
% % Convert volume to surface with vol2sur
% reg_file = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/fsaverage/mri.2mm/reg.2mm.dat';
% unix(['mri_vol2surf --src ' filename '_size_FDR_001_raw_statmap.nii --out ./' filename '_size_FDR_001_raw_statmap_surf_lh.nii --reg ' reg_file ' --trgsubject fsaverage --hemi lh' ])
% unix(['mri_vol2surf --src ' filename '_size_FDR_001_thresh_statmap.nii --out ./' filename '_size_FDR_001_thresh_statmap_surf_lh.nii --reg ' reg_file ' --trgsubject fsaverage --hemi lh' ])
% unix(['mri_vol2surf --src ' filename '_size_FDR_001_pmap.nii --out ./' filename '_size_FDR_001_pmap_surf_lh.nii --reg ' reg_file ' --trgsubject fsaverage --hemi lh' ])
% 
% unix(['mri_vol2surf --src ' filename '_size_FDR_001_raw_statmap.nii --out ./' filename '_size_FDR_001_raw_statmap_surf_rh.nii --reg ' reg_file ' --trgsubject fsaverage --hemi rh' ])
% unix(['mri_vol2surf --src ' filename '_size_FDR_001_thresh_statmap.nii --out ./' filename '_size_FDR_001_thresh_statmap_surf_rh.nii --reg ' reg_file ' --trgsubject fsaverage --hemi rh' ])
% unix(['mri_vol2surf --src ' filename '_size_FDR_001_pmap.nii --out ./' filename '_size_FDR_001_pmap_surf_rh.nii --reg ' reg_file ' --trgsubject fsaverage --hemi rh' ])
% 

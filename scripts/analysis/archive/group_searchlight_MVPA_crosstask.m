%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to perform crosstask searchlight MVPA analysis on the
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

conds_perblock_str = 'blocktrials'; % blocktrials or 16trials
conds_perblock = 1; %16;

cond_train_str = 'WM-SMC';
cond_train = [1,2,3,4,5,6];
cond_train_labels = [1,1,0,1,1,0];
cond_test_str = 'SMC-F';
cond_test = [3,6,10];
cond_test_labels = [1,1,0];


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

if any(ismember([cond_train, cond_test], 10)) && any(ismember(subjCodes, {'MM', 'PP'}))
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
N_runs = 4; % all subjs have 4 runs (except PP has 6, this will be accounted for)

N_trials_perrun_train = length(cond_train); 
N_trials_persubj_train = N_trials_perrun_train * N_runs; 
x_train = nan((N_subjs * N_trials_persubj_train)+(trials_mod*N_trials_perrun_train), sum(cortex_mask_reshaped,'all'));
y_train = nan((N_subjs * N_trials_persubj_train)+(trials_mod*N_trials_perrun_train), 1);
subjlist_train = nan((N_subjs * N_trials_persubj_train)+(trials_mod*N_trials_perrun_train), 1);

N_trials_perrun_test = length(cond_test); 
N_trials_persubj_test = N_trials_perrun_test * N_runs;
x_test = nan((N_subjs * N_trials_persubj_test)+(trials_mod*N_trials_persubj_test), sum(cortex_mask_reshaped,'all'));
y_test = nan((N_subjs * N_trials_persubj_test)+(trials_mod*N_trials_persubj_test), 1);
subjlist_test = nan((N_subjs * N_trials_persubj_test)+(trials_mod*N_trials_perrun_test), 1);

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

        % Get coefs for relevant conditions
        coefs_train = coefs.vol(:,:,:,cond_train);
        coefs_test = coefs.vol(:,:,:,cond_test);
        coefs_train_reshaped = reshape(coefs_train, size(coefs_train,1)*size(coefs_train,2)*size(coefs_train,3), length(cond_train));
        coefs_test_reshaped = reshape(coefs_test, size(coefs_test,1)*size(coefs_test,2)*size(coefs_test,3), length(cond_test));
        coefs_train_inmask = coefs_train_reshaped(cortex_mask_reshaped==1,:);
        coefs_test_inmask = coefs_test_reshaped(cortex_mask_reshaped==1,:);

        % Z-score coefs
        coefs_train_zscored = zscore(coefs_train_inmask'); % This is z-scored each column (voxel) separately
        coefs_test_zscored = zscore(coefs_test_inmask'); % This is z-scored each column (voxel) separately

        % Store data
        index_start_train = ((count-1)*N_trials_perrun_train + 1);
        indices_train = index_start_train:(index_start_train+N_trials_perrun_train-1);
        x_train(indices_train,:) = coefs_train_zscored;
        y_train(indices_train) = cond_train_labels; % always the same because coefs are stored in the same condition order (1 = WM, 0 = SMC)
        subjlist_train(indices_train) = repmat(ss, N_trials_perrun_train,1);

        index_start_test = ((count-1)*N_trials_perrun_test + 1);
        indices_test = index_start_test:(index_start_test+N_trials_perrun_test-1);
        x_test(indices_test,:) = coefs_test_zscored;
        y_test(indices_test) = cond_test_labels; % always the same because coefs are stored in the same condition order (1 = WM, 0 = SMC)
        subjlist_test(indices_test) = repmat(ss, N_trials_perrun_test,1);

    end
    disp(['Finished subj ' subjCode])
end

%save('trial_coefs_WM_SMC_pertrial.mat', 'x', 'y', 'subjlist', 'runlist');

%load('/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/MVPA_results/trial_coefs_SMC_F.mat', 'x', 'y', 'subjlist', 'runlist');

%% Run searchlight analysis
CVfold_train = subjlist_train; % information for leave-one-subject-out
CVfold_test = subjlist_test;
dist = 2; % radius of searchlight, in voxels. radius of 2 gives a small sphere that is 33 voxels
metric = 'Pearson'; % will use correlation between choice and prediction to gauge performance
method = 'PLS1'; 

% Remove nans (this should only happen if a full run is bad)
if any(isnan(x_train),'all')
    nanmask_train = any(isnan(x_train),2);
    check_nans_train = sum(nanmask_train);
    assert(mod(check_nans_train, N_trials_perrun_train)==0, 'sum of nans are not a multiple of trials per run');
    x_train = x_train(~nanmask_train,:);
    y_train = y_train(~nanmask_train);
    CVfold_train = CVfold_train(~nanmask_train);
end
if any(isnan(x_test),'all')
    nanmask_test = any(isnan(x_test),2);
    check_nans_test = sum(nanmask_test);
    assert(mod(check_nans_test, N_trials_perrun_test)==0, 'sum of nans are not a multiple of trials per run');
    x_test = x_test(~nanmask_test,:);
    y_test = y_test (~nanmask_test);
    CVfold_test = CVfold_test(~nanmask_test);
end

tic;
measure = Afmri.searchlight_crosstask(x_train, x_test, y_train, y_test, CVfold_train, CVfold_test, cortex_mask, dist, metric, method); % nested leave-one-subject-out cross-validation to find searchlights that can predict Y
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
MRIwrite(measure_data, [results_dir 'crosstask_corr_measure_' cond_train_str '_' cond_test_str '_' conds_perblock_str '.nii']);

mask_data = template;
mask_data.vol(cortex_mask) = 1;
mask_data.vol(~cortex_mask) = 0;
MRIwrite(mask_data, [results_dir 'crosstask_mask_' cond_train_str '_' cond_test_str '_' conds_perblock_str '.nii']);

%% Run TFCE permutation testing
output_file = [results_dir 'crosstask_' cond_train_str '_' cond_test_str '_' conds_perblock_str '_MVPA_FWE'];
input_data = [results_dir 'crosstask_corr_measure_' cond_train_str '_' cond_test_str '_' conds_perblock_str '.nii'];
mask = [results_dir 'crosstask_mask_' cond_train_str '_' cond_test_str '_' conds_perblock_str '.nii'];
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

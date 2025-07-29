%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to perform searchlight MVPA analysis on the
% trifloc localizer data to predict presence of each condition 
% Tom Possidente - July 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Notes:
% - should start with individual-level predictions
% - may have to split blocks into 2 or more to get enough samples per
% subj/run

addpath(genpath('/projectnb/somerslab/tom/functions/'));
addpath(genpath('/projectnb/somerslab/tom/ArthurfMRI-main/'));
ccc;

%% Set up key variables and paths
% Load subj codes
experiment_name = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
% subjCodes{1} = subjDf_cut.subjCode(~ismember(subjDf_cut.subjCode, {'AH', 'SL', 'RR'})); % N=21 (2 are missing fixation condition)
% subjCodes{2} = subjDf_cut.subjCode(~ismember(subjDf_cut.subjCode, {'AH', 'SL', 'RR', 'PP', 'MM'})); % N=19 (all have fixation condition)
% N_subjs{1} = length(subjCodes{1});
% N_subjs{2} = length(subjCodes{2});
subjCodes = {'GG'};
N_subjs = length(subjCodes);

basedir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';

conditions = 1:4; % all 4 visual active sections  
% localizer condition order (sections_per_block=2):
% 1-4 = vA (block 1, section 1; block 1, section 2; block 2, section 1; block 2, section 2)
% 5-6 = vP (block 1, section 1; block 1, section 2)
% 7-10 = aA (block 1, section 1; block 1, section 2; block 2, section 1; block 2, section 2)
% 11-12 = aP (block 1, section 1; block 1, section 2)
% 13-16 = tA (block 1, section 1; block 1, section 2; block 2, section 1; block 2, section 2)
% 17-18 = tP (block 1, section 1; block 1, section 2)
% 19-20 = f (block 1, section 1; block 1, section 2)

%% For each subj/run, load GLM results and perform MVPA searchlight

for ss = 1:N_subjs
    subjCode = subjCodes{ss};

    % Get # of runs for subj
    files = {dir([basedir subjCode '/localizer/']).name};
    N_runs = sum(contains(files, '00'));
    
    for rr = 1:N_runs
        
        % Load GLM coefficients
        coefs = MRIread([basedir subjCode '/localizer/localizer_block_regression_0sm_run' num2str(rr) '/beta.nii.gz']);

        % Load param file to get coef order
        params = readmatrix([basedir subjCode '/localizer/00' num2str(rr) '/localizer_condition_timing_perblock.para'], 'FileType','text');

        % Get coefs for relevant conditions
        cond_inds = find(ismember(params(:,2), conditions));
        coefs = coefs.vol(:,:,:,cond_inds);
        coefs_reshape = reshape(coefs, size(coefs,1)*size(coefs,2)*size(coefs,3), length(conditions));

        % 


    end

end








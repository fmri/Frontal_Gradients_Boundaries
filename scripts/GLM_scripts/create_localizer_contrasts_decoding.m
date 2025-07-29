%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to run block-level regressions for
%%% MVPA decoding purposes on the trifloc localizer data
%%% Tom Possidente - July 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/functions/');
ccc;

%% Set up directories and subj info
experiment_name = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
% subjCodes{1} = subjDf_cut.subjCode(~ismember(subjDf_cut.subjCode, {'AH', 'SL', 'RR'})); % N=21 (2 are missing fixation condition)
% subjCodes{2} = subjDf_cut.subjCode(~ismember(subjDf_cut.subjCode, {'AH', 'SL', 'RR', 'PP', 'MM'})); % N=19 (all have fixation condition)
% N_subjs{1} = length(subjCodes{1});
% N_subjs{2} = length(subjCodes{2});
subjCodes = {'GG'};
N_subjs = 1;

data_dirbase = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';

analysis_name = 'localizer_block_regression_0sm_run1';

TR = 2;
para_name = 'localizer_condition_timing_perblock.para';

funcstem = 'fmcpr_tu.siemens.sm0.mni305.2mm.nii.gz'; % registered, fieldmaps applied, motion corrected, slice time corrected volumes (no smoothing)
rlf_name = 'localizer_contrasts_runlistfile_run1.txt';

sections_per_block = 2;
refeventdur = 32/sections_per_block;
nconditions = 20; 

run_mkanalysis = true;
run_mkcontrast = true;

% localizer condition orders (1way)
% 1 = vA
% 2 = vP
% 3 = aA
% 4 = aP
% 5 = tA
% 6 = tP
% 7 = f


%% Create analysis/contrast info

% using -event-related (used for block or event-related designs)
% Using polyfit 2, spmhrf 1 to use 1 derivative of the hrf.
% mcextreg option will use 12 motion regressors (3 translation 3 rotation
% and derivatives)

if run_mkanalysis
    unix(['mkanalysis-sess -a ' analysis_name ' -funcstem ' funcstem ... 
        ' -mni305 2 -fsd localizer -event-related -paradigm ' para_name ...  
        ' -nconditions ' num2str(nconditions) ' -refeventdur ' num2str(refeventdur) ' -TR ' num2str(TR) ... % refeventdur should be number of seconds for a block (or whatever unit of time each condition is)
        ' -polyfit 2 -spmhrf 1 -mcextreg -runlistfile ' rlf_name ' -per-run -force']) 
end


%% run make contrast

% localizer condition order (sections_per_block=2):
% 1-4 = vA (block 1, section 1; block 1, section 2; block 2, section 1; block 2, section 2)
% 5-6 = vP (block 1, section 1; block 1, section 2)
% 7-10 = aA (block 1, section 1; block 1, section 2; block 2, section 1; block 2, section 2)
% 11-12 = aP (block 1, section 1; block 1, section 2)
% 13-16 = tA (block 1, section 1; block 1, section 2; block 2, section 1; block 2, section 2)
% 17-18 = tP (block 1, section 1; block 1, section 2)
% 19-20 = f (block 1, section 1; block 1, section 2)

if run_mkcontrast
    for cc = 1:4
        unix(['mkcontrast-sess -analysis ' analysis_name ' -contrast ' ...
            'vA' num2str(cc) ' -a ' num2str(cc)])
    end
end

%% Loop through subjs
for ss = 1:N_subjs

    subjCode = subjCodes{ss};

    % run glm
    unix(['selxavg3-sess -s ' subjCode ' -d ' data_dirbase ' -analysis ' analysis_name...
        ' -no-preproc -overwrite'])

end

















%%
%%%%%%%%%% Run for subjs RR MM PP (no fixation condition) %%%%%%%%%%
%%

analysis_name_lh_nofix = 'localizer_contrasts_0sm_nofix_lh';
analysis_name_rh_nofix = 'localizer_contrasts_0sm_nofix_rh';
rlf_name = 'localizer_contrasts_runlistfile.txt';
subjCodes = {'MM', 'PP'};

run_mkanalysis = true;
run_mkcontrast = true;

if run_mkanalysis
    unix(['mkanalysis-sess -a ' analysis_name_lh_nofix ' -funcstem ' funcstem_lh ...
        ' -surface fsaverage lh -fsd localizer -event-related -paradigm ' para_name ...
        ' -nconditions 6 -refeventdur 32 -TR ' num2str(TR) ...
        ' -polyfit 1 -spmhrf 0 -mcextreg -runlistfile ' rlf_name ' -per-run -force'])

    unix(['mkanalysis-sess -a ' analysis_name_rh_nofix ' -funcstem ' funcstem_rh ...
        ' -surface fsaverage rh -fsd localizer -event-related -paradigm ' para_name ...
        ' -nconditions 6 -refeventdur 32 -TR ' num2str(TR) ...
        ' -polyfit 1 -spmhrf 0 -mcextreg -runlistfile ' rlf_name ' -per-run -force'])
end


if run_mkcontrast

        unix(['mkcontrast-sess -analysis ' analysis_name_lh_nofix ' -contrast ' ...
            'vA-aA -a 1 -c 3'])
        unix(['mkcontrast-sess -analysis ' analysis_name_rh_nofix ' -contrast ' ...
            'vA-aA -a 1 -c 3'])

        unix(['mkcontrast-sess -analysis ' analysis_name_lh_nofix ' -contrast ' ...
            'vP-aP -a 2 -c 4'])
        unix(['mkcontrast-sess -analysis ' analysis_name_rh_nofix ' -contrast ' ...
            'vP-aP -a 2 -c 4'])

        unix(['mkcontrast-sess -analysis ' analysis_name_lh_nofix ' -contrast ' ...
            'vA-vP -a 1 -c 2'])
        unix(['mkcontrast-sess -analysis ' analysis_name_rh_nofix ' -contrast ' ...
            'vA-vP -a 1 -c 2'])

        unix(['mkcontrast-sess -analysis ' analysis_name_lh_nofix ' -contrast ' ...
            'aA-aP -a 3 -c 4'])
        unix(['mkcontrast-sess -analysis ' analysis_name_rh_nofix ' -contrast ' ...
            'aA-aP -a 3 -c 4'])

        unix(['mkcontrast-sess -analysis ' analysis_name_lh_nofix ' -contrast ' ...
            'V-A -a 1 -a 2 -c 3 -c 4'])
        unix(['mkcontrast-sess -analysis ' analysis_name_rh_nofix ' -contrast ' ...
            'V-A -a 1 -a 2 -c 3 -c 4'])

        unix(['mkcontrast-sess -analysis ' analysis_name_lh_nofix ' -contrast ' ...
            'vAaA-vPaP -a 1 -a 3 -c 2 -c 4'])
        unix(['mkcontrast-sess -analysis ' analysis_name_rh_nofix ' -contrast ' ...
            'vAaA-vPaP -a 1 -a 3 -c 2 -c 4'])

end

for ss = 1:length(subjCodes)

    subjCode = subjCodes{ss};

    % run glm
    unix(['selxavg3-sess -s ' subjCode ' -d ' data_dir ' -analysis ' analysis_name_lh_nofix ...
        ' -no-preproc -overwrite'])
    unix(['selxavg3-sess -s ' subjCode ' -d ' data_dir ' -analysis ' analysis_name_rh_nofix ...
        ' -no-preproc -overwrite'])

end

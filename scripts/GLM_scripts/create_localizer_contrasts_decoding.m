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
subjCodes = subjDf_cut.subjCode(~ismember(subjDf_cut.subjCode, {'AH', 'SL', 'RR', 'PP', 'MM'})); % N=19 (all have fixation condition)
subjCodes = {'LA'}
N_subjs = length(subjCodes);

data_dirbase = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';

analysis_name = 'localizer_block_regression_0sm_run1_16sections';

TR = 2;
para_name = 'localizer_condition_timing_perblock_16sections.para';

funcstem = 'fmcpr_tu.siemens.sm0.mni305.2mm.nii.gz'; % registered, fieldmaps applied, motion corrected, slice time corrected volumes (no smoothing)
rlf_name = 'localizer_contrasts_runlistfile_run1.txt';

sections_per_block = 16;
refeventdur = 32/sections_per_block;
nconditions = 10*sections_per_block; 

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
% Using polyfit 2, spmhrf 0 to use 0 derivative of the hrf.
% -nuisreg mcprextreg 6 will use 6 motion regressors

if run_mkanalysis
    unix(['mkanalysis-sess -a ' analysis_name ' -funcstem ' funcstem ... 
        ' -mni305 2 -fsd localizer -event-related -paradigm ' para_name ...  
        ' -nconditions ' num2str(nconditions) ' -refeventdur ' num2str(refeventdur) ' -TR ' num2str(TR) ... % refeventdur should be number of seconds for a block (or whatever unit of time each condition is)
        ' -polyfit 2 -spmhrf 0 -nuisreg mcprextreg 6 -runlistfile ' rlf_name ' -per-run -force']) 
end


%% run make contrast

% localizer condition order (sections_per_block=16):
% 1-32 = vA 
% 33-48 = vP 
% 49-80 = aA 
% 81-96 = aP 
% 97-128 = tA 
% 129-144 = tP 
% 145-160 = f 
% localizer condition order (sections_per_block=4):
% 1-8 = vA (2 blocks, 4 sections each)
% 9-12 = vP (1 block, 4 sections each)
% 13-20 = aA (2 blocks, 4 sections each)
% 21-24 = aP (1 block, 4 sections each)
% 25-32 = tA (2 blocks, 4 sections each)
% 33-36 = tP (1 block, 4 sections each)
% 37-40 = f (1 block, 4 sections each)

condition_names = [arrayfun(@(x) 'vA', 1:sections_per_block*2, 'UniformOutput',false),...
                   arrayfun(@(x) 'vP', 1:sections_per_block, 'UniformOutput',false),...
                   arrayfun(@(x) 'aA', 1:sections_per_block*2, 'UniformOutput',false),...
                   arrayfun(@(x) 'aP', 1:sections_per_block, 'UniformOutput',false),...
                   arrayfun(@(x) 'tA', 1:sections_per_block*2, 'UniformOutput',false),...
                   arrayfun(@(x) 'tP', 1:sections_per_block, 'UniformOutput',false),...
                   arrayfun(@(x) 'f', 1:sections_per_block, 'UniformOutput',false)]';
cond_inds = [1:96, 145:160];

if run_mkcontrast
    parfor cc = 1:sections_per_block*7
        cond_ind = cond_inds(cc);
        unix(['mkcontrast-sess -analysis ' analysis_name ' -contrast ' ...
            condition_names{cond_ind} num2str(cond_ind) ' -a ' num2str(cond_ind)])
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

analysis_name_nofix = 'localizer_block_regression_0sm_run6_16section_nofix';
rlf_name = 'localizer_contrasts_runlistfile_run6.txt';
nconditions = 12*sections_per_block; 
subjCodes = {'PP'};

run_mkanalysis = true;
run_mkcontrast = true;

if run_mkanalysis
    unix(['mkanalysis-sess -a ' analysis_name_nofix ' -funcstem ' funcstem ... 
        ' -mni305 2 -fsd localizer -event-related -paradigm ' para_name ...  
        ' -nconditions ' num2str(nconditions) ' -refeventdur ' num2str(refeventdur) ' -TR ' num2str(TR) ... % refeventdur should be number of seconds for a block (or whatever unit of time each condition is)
        ' -polyfit 2 -spmhrf 0 -nuisreg mcprextreg 6 -runlistfile ' rlf_name ' -per-run -force']) 
end


if run_mkcontrast
    condition_names = [arrayfun(@(x) 'vA', 1:sections_per_block*2, 'UniformOutput',false),...
                   arrayfun(@(x) 'vP', 1:sections_per_block*2, 'UniformOutput',false),...
                   arrayfun(@(x) 'aA', 1:sections_per_block*2, 'UniformOutput',false),...
                   arrayfun(@(x) 'aP', 1:sections_per_block*2, 'UniformOutput',false),...
                   arrayfun(@(x) 'tA', 1:sections_per_block*2, 'UniformOutput',false),...
                   arrayfun(@(x) 'tP', 1:sections_per_block*2, 'UniformOutput',false)]';
    parfor cc = 1:sections_per_block*8
        unix(['mkcontrast-sess -analysis ' analysis_name_nofix ' -contrast ' ...
            condition_names{cc} num2str(cc) ' -a ' num2str(cc)])
    end
end

parfor ss = 1:length(subjCodes)

    subjCode = subjCodes{ss};

    % run glm
    unix(['selxavg3-sess -s ' subjCode ' -d ' data_dirbase ' -analysis ' analysis_name_nofix...
        ' -no-preproc -overwrite'])

end

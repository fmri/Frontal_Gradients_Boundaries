%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to create runlist txt files for odd and
%%% even runs for the GLM to use to pick which runs to include in analysis
%%% during cross-validation
%%% 
%%% Tom Possidente - September 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/');
ccc;

%% Set up directories and subj info

projectDir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/';

subjCodes = {'MK','AB', 'AD','LA','AE','TP','NM','AF','AG','AI','GG','UV','PQ','KQ','LN','RT','PT','PL','NS','SL'};
N = length(subjCodes);
odd_run_nums = {'001'; '003'};
even_run_nums = {'002'; '004'};

%% Loop over subjs and create files
for nn = 1:N
    subjCode = subjCodes{nn};
    filepath_base = [projectDir 'data/unpacked_data_nii_fs_localizer/' subjCode '/localizer/'];
    writecell(odd_run_nums, [filepath_base 'localizer_contrasts_runlistfile_odds.txt']);
    writecell(even_run_nums, [filepath_base 'localizer_contrasts_runlistfile_evens.txt']);
end

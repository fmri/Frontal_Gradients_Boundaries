%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to take preprocessed subject space
%%% volumes and transform them to MNI305 space
%%% Tom Possidente - July 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/functions/');
ccc;

%% Set key variables and paths
experiment_name = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
N_subjs = length(subjCodes);

basedir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer';
funcstem = 'fmcpr_tu.siemens';

%% Loop through subjs and transform volumes to mni305
parfor ss = 1:N_subjs

    subjCode = subjCodes{ss};

    unix(['rawfunc2tal-sess -i fmcpr_tu.siemens -fwhm 0 -s ' subjCode ' -d ' basedir ' -fsd localizer -per-run -no-subcort-mask'])

end













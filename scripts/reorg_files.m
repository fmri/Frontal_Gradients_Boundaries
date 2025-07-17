
addpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/');
ccc;

experiment_name = 'spacetime';

projectDir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/';

subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjectsDir_src = [projectDir, 'data/unpacked_data_nii_fs_localizer/'];
hemis = {'lh', 'rh'};

%% Loop over subjs
for ss = 1:length(subjCodes)
    subjCode = subjCodes{ss};
    subjDir = [subjectsDir_src subjCode '/localizer/'];

    for hh = 1:length(hemis)
        hemi = hemis{hh};
        src = [subjDir 'localizer_contrasts_extra_nofix_' hemi '/vA-aA/'];

        unix(['cp -r ' src ' ' subjDir 'localizer_contrasts_' hemi '/']);
    end
end












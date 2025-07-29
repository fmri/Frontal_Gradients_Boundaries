%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to perform a group level activity
%%% analysis on the spacetime localizer data sensory drive and WM contrasts
%%% for visual and auditory stimulation/tasks.
%%%
%%% Tom Possidente - December 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% load subjIDs
experiment_name = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
%subjCodes = subjDf_cut.subjCode(~ismember(subjDf_cut.subjCode, {'AH', 'SL', 'RR'}));
subjCodes{1} = subjDf_cut.subjCode(~ismember(subjDf_cut.subjCode, {'AH', 'SL', 'RR'})); % active contrasts
subjCodes{2} = subjDf_cut.subjCode(~ismember(subjDf_cut.subjCode, {'AH', 'SL', 'RR', 'PP', 'MM'})); % passive contrasts
N_subjs{1} = length(subjCodes{1});
N_subjs{2} = length(subjCodes{2});

fs_num = 163842;

%% Loop through subjs and contrasts and construct mri_concat command
contrasts = {'f-vP', 'f-aP', 'aAaP-f', 'vAvP-f', 'aA-aP', 'vA-vP', 'vAaA-vPaP', 'aPvP-f', 'V-A'};
%contrasts = {'aA-aP', 'vA-vP', 'tA-tP'};
% contrasts = {'vAaA-vPaP', 'aPvP-f'};
%contrasts = {'V-A'};
subjCodes_list = [2,2,2,2,1,1,1,2,1];
reverse_contrast = [true, true, false, false, false, false, false, false, false];

design = 'design.con'; % design.con for active contrasts, design_neg.con for passive contrasts
N_contrasts = length(contrasts);
subj_path = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';
path_suffix = '/localizer/localizer_contrasts_';
hemis = {'lh', 'rh'};

% already ran this
for hh = 1:length(hemis)
    hemi = hemis{hh};
    for cc = 1:N_contrasts
        cmd = 'mri_concat ';
        contrast = contrasts{cc};
        for ss = 1:N_subjs{subjCodes_list(cc)}
            subjCode = subjCodes{subjCodes_list(cc)}{ss};
            path = [subj_path subjCode path_suffix hemi '/' contrast '/ces.nii.gz'];
            cmd = [cmd path ' '];
        end
        cmd = [cmd '--o ' hemi '.ces.localizer_groupavg_' contrast '.nii'];
        unix(cmd);
    end
end

%% Loop through concatenated data and run glm for each contrast and hemisphere
groupavg_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/grouplevel/';

% % Already ran this
for cc = 1:N_contrasts
    contrast = contrasts{cc};
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        outdir = ['/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/grouplevel/' hemi '.ces.localizer_groupavg_' contrast '.glmres'];
        unix(['mri_glmfit --y ' groupavg_dir hemi '.ces.localizer_groupavg_' contrast '.nii --surf fsaverage ' hemi ' --osgm --glmdir ' outdir ' --cortex']);
    end
end

%% Multiple comparison correction (TFCE)
% FSL's randomise has a TFCE option but does not support surface data, but
% PALM does, so using that 

addpath(genpath('/projectnb/somerslab/tom/palm-alpha119'));
res_path = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/group_contrast_TFCE/';
fsdir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/fsaverage/surf/';
tfce_iters = '1000';
dmats = {'design21', 'design19'};

masks = {'/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/lh.frontal_cortex_spmreshaped.nii',
         '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/rh.frontal_cortex_spmreshaped.nii'};

parfor cc = 1:N_contrasts
    contrast = contrasts{cc};
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        input = [res_path hemi '.ces.localizer_groupavg_' contrast '.nii'];
        outdir = [res_path contrast '/' hemi '/'];
        if reverse_contrast(cc)
            splt_contrast = split(contrast,'-');
            contrast_str = [splt_contrast{2} '-' splt_contrast{1}];
            outdir = [res_path contrast_str '/' hemi '/'];
        end
        if ~isfolder(outdir)
            mkdir(outdir)
        end
        
        % Reshape input data because FSL doesn't work nicely with freesurfer surface outputs
        spm_nii_path = [res_path hemi '.ces.localizer_groupavg_' contrast '_reshaped.nii'];
        fs2spm_nii(input, spm_nii_path, reverse_contrast(cc));
        
        % Run tfce with PALM
        palm('-i', spm_nii_path, '-d', [res_path dmats{subjCodes_list(cc)} '.mat'], '-t' ,[res_path design], '-s',...
            [fsdir hemi '.pial'], '-m', masks{hh}, '-n', tfce_iters, '-ise', '-T', '-tfce2D', '-logp', '-o', outdir)

    end
end

%% PALM is outputting in nifti-2 format for some reason, so use nifti22nifti1.py to convert back into nifti-1 so matlab/freesurfer can use the files
for cc = 1:N_contrasts
    contrast = contrasts{cc};
    if reverse_contrast(cc)
        splt_contrast = split(contrast,'-');
        contrast = [splt_contrast{2} '-' splt_contrast{1}];
    end
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        input = [res_path contrast '/' hemi '/' '_tfce_tstat_fwep.nii'];
        output = [res_path contrast '/' hemi '/' '_tfce_tstat_fwep_nifti1.nii'];
        function_path = '/projectnb/somerslab/tom/functions/nifti22nifti1.py';
        cmd = sprintf('python3 %s %s %s', function_path, input, output);
        unix(['module load miniconda && conda activate tom_env && ' cmd])
    end
end



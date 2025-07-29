%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to call freeview_screenshots to take
% screenshots of the lefco normalized difference between contrasts for
% individual subjects
%
% Tom Possidente August 1 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/'))
ccc;

%% Get subjIDs
experiment_name = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjCodes = subjDf_cut.subjCode(~ismember(subjDf_cut.subjCode, {'AH', 'SL', 'RR', 'PP', 'MM'}));
N = length(subjCodes);


%% Set ROI and contrast lists

contrast = 'vAaA-vPaP-aPvP-f';

lh_ROI_paths = {};
rh_ROI_paths = {};
contrasts = cell(N,1);
stat_file = cell(N,2);

for ss = 1:N % loop through subjs
    subjCode = subjCodes{ss};
    
    contrasts{ss,1} = contrast;
    
    stat_file{ss,1} = [subjCode '_lh_' contrast '_norm_frontal.nii.gz'];
    stat_file{ss,2} = [subjCode '_rh_' contrast '_norm_frontal.nii.gz'];
    % lh_ROI_paths{ss,1} = '';
    % lh_ROI_paths{ss,2} = '';
    % 
    % rh_ROI_paths{ss,1} = '';
    % rh_ROI_paths{ss,2} = '';

end

%% Actually run freeview_screenshots()
ctable = [];
lh_ROI_paths = {};
rh_ROI_paths = {};
savedir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/figures_images/norm_gradients_indiv/';
func_path = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/';
func_folder = 'data';
analysis_name = 'gradient_niftis';
freeview_screenshots(subjCodes, contrasts, lh_ROI_paths, rh_ROI_paths, ctable, true, savedir, func_path, func_folder,...
    [], analysis_name, stat_file,'0.001,1', [], [], [], true);
crop_ppt_fv_images(subjCodes, contrasts)

% subjIDs, contrast_list, lh_label_list, rh_label_list, colortable, use_fsaverage, save_dir, func_path, func_folder, recon_dir, 
% analysis_name, stat_file, overlayThreshold, label_opacity, ss_suffix, ss_on

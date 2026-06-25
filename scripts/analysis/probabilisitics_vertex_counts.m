%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to calculate the number of vertices in
%%% the SMC and WM versions of the probabilistic maps
%%%
%%% Tom Possidente - June 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables
ROI_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/probabilistic_allROIs/';
unpermuted_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/SMC_WM_nonpermuted_probabilistics/';

hemis = {'lh', 'rh'};
contrasts = {'WM', 'SMC'};
modalities = {'visual', 'auditory'};
thresh = 5;

ROI_names = {{'aMFG', 'midIFS', 'aINS', 'preSMA', 'inf_lat_frontal', 'sup_lat_frontal', ...
    'aIPS', 'VOT', 'cIPS', 'LOT', 'pIPS', 'VO', 'DO'},...
    {'aMFG', 'midIFS', 'aINS', 'preSMA', 'inf_lat_frontal', 'sup_lat_frontal', ...
    'parietal_opercular', 'aIPS', 'ms_post_STSG', 'cIPS'} }; % ROI names for visual then auditory modality

ROI_names = {{'midIFS', 'aINS', 'preSMA', 'inf_lat_frontal', 'sup_lat_frontal', ...
    'aIPS'},...
    {'midIFS', 'aINS', 'preSMA', 'inf_lat_frontal', 'sup_lat_frontal', ...
     'aIPS'} }; % ROI names for visual then auditory modality

%% Loop through probabilistics and count vertices
n_vertices_vis = nan(2,2,length(ROI_names{1})); % hemisphere, contrast, ROI
n_vertices_aud = nan(2,2,length(ROI_names{2}));
n_vertices = {n_vertices_vis, n_vertices_aud};
n_verts = table();
perc_verts = table();
vert_inds = 1:163842;

for hh = 1:2
    hemi = hemis{hh};
    for mm = 1:2
        modality = modalities{mm};
        for cc = 1:2
            contrast = contrasts{cc};
            probmap = MRIread([unpermuted_dir hemi '_original_' contrast  '_' modality '.nii']);
            ROI_inds_good = probmap.vol >= thresh;
            for rr = 1:length(ROI_names{mm})
                ROI_name = ROI_names{mm}{rr};
                ROI = readtable([ROI_dir hemi '.' ROI_name '_prob_thresh5.label'], 'FileType','text');
                ROI_inds = ROI{:,1};
                ROI_probs_threshed = probmap.vol(ROI_inds_good & ismember(vert_inds,ROI_inds+1));
                n_verts(end+1,:) = {ROI_name, hemi, modality, contrast, length(ROI_probs_threshed)};
                perc_verts(end+1,:) = { ROI_name, hemi, modality, contrast, length(ROI_probs_threshed)/length(ROI_inds) }
            end
        end
    end
end

n_verts.Properties.VariableNames = {'ROI', 'hemi', 'modality', 'contrast', 'num_vertices'}
perc_verts.Properties.VariableNames = {'ROI', 'hemi', 'modality', 'contrast', 'perc_verts'}
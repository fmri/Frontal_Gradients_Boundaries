%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to calculate dice coefficients between 
%%% sensory and WM probabilistics and frontal domain general MMP1 parcels
%%%
%%% Tom Possidente - June 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables
hemis = {'lh', 'rh'};
N_hemis = length(hemis);

modalities = {'auditory', 'visual'};
N_modalities = length(modalities);

contrasts = {'SMC', 'WM'};
N_contrasts = length(contrasts);

MD_ROI_names = {'8BM_ROI', '8C_ROI', 'IFJp_ROI', 'p9-46v_ROI', 'a9-46v_ROI', 'i6-8_ROI', 'AVI_ROI'};
N_MD_ROIs = length(MD_ROI_names);

prob_ROI_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/probabilistic_allROIs/';
prob_ROIs = {'preSMA', 'inf_lat_frontal', 'aINS', 'sup_lat_frontal', 'midIFS'};
N_prob_ROIs = length(prob_ROIs);

N_vertices = 163842;

%% Extract MD labels from annotation files
MD_annot_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/fsaverage/label/';
MD_ROI_masks = nan(N_hemis, N_MD_ROIs, N_vertices);

for hh = 1:N_hemis
    [vertices, label, colortable] = read_annotation([MD_annot_dir hemis{hh} '.HCP-MMP1.annot']);
    for rr = 1:N_MD_ROIs
        label_mask = strcmpi(colortable.struct_names, [hemis{hh}(1) '_' MD_ROI_names{rr}]);
        assert(sum(label_mask)==1, '0 or >1 label found with same name in annotation file, should not happen');
        label_num = colortable.table(label_mask,5);
        MD_ROI_masks(hh,rr,:) = label == label_num;
    end
end

MD_ROI_masks = squeeze(sum(MD_ROI_masks,2));

%% Loop over hemis/ROIs/modalities/contrasts and calculate dice coefs
probmap_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/SMC_WM_nonpermuted_probabilistics/';
dice_coefs = nan(N_hemis, N_modalities, N_contrasts, N_prob_ROIs);

for cc = 1:N_contrasts
    contrast = contrasts{cc};
    for mm = 1:N_modalities
        modality = modalities{mm};
        for hh = 1:N_hemis
            hemi = hemis{hh};
            probmap = MRIread([probmap_dir hemi '_original_' contrast '_' modality '.nii']);
            probmap_mask = probmap.vol>=5;
            MD_ROI_mask = MD_ROI_masks(hh,:);
            for rr = 1:N_prob_ROIs
                prob_ROI = readtable([prob_ROI_dir hemi '.' prob_ROIs{rr} '_prob_thresh5.label'], 'FileType','text');
                prob_ROI_mask = ismember(1:N_vertices, (prob_ROI{:,1}+1));
                probmask_inROI = prob_ROI_mask & probmap_mask;
                MD_ROI_inROI = MD_ROI_mask & prob_ROI_mask;
                overlap = probmask_inROI & MD_ROI_mask;
                dice_coefs(hh,mm,cc,rr) = 2*sum(overlap) / (sum(probmask_inROI) + sum(MD_ROI_inROI));
                disp([prob_ROIs{rr} ' ' modality ' ' contrast ' ' hemi ': ' num2str(dice_coefs(hh,mm,cc,rr))]);
            end
        end
    end
end


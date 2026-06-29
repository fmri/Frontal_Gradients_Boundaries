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

dice_coefs = nan(N_hemis, N_modalities, N_contrasts);

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
dice_coefs = nan(N_hemis, N_modalities, N_contrasts);

for cc = 1:N_contrasts
    contrast = contrasts{cc};
    for mm = 1:N_modalities
        modality = modalities{mm};
        for hh = 1:N_hemis
            hemi = hemis{hh};
            probmap = MRIread([probmap_dir hemi '_original_' contrast '_' modality '.nii']);
            probmap_mask = probmap.vol>=5;
            probmap_vertices = sum(probmap_mask);
            MD_ROI_mask = MD_ROI_masks(hh,:);
            overlap = probmap_mask & MD_ROI_mask;
            dice_coefs(hh,mm,cc) = 2*sum(overlap) / (sum(probmap_vertices) + sum(MD_ROI_mask));
            disp([modality ' ' contrast ' ' hemi ': ' num2str(dice_coefs(hh,mm,cc))]);
        end
    end
end



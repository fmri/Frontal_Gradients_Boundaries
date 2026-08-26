%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to find and select the Glassar atlas
%%% parcels that overlay with the probabilistic activity for each contrast
%%% of interest so that we can display them in freeview to make figures.
%%%
%%% Tom Possidente - May 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables
label_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/probabilistic_allROIs/';
ROI_names = {'midIFS', 'aINS', 'preSMA', 'inf_lat_frontal', 'sup_lat_frontal', 'aIPS'};
[vertices, label_codes{1}, colortables{1}] = read_annotation('/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/fsaverage/label/lh.HCP-MMP1.annot');
[~, label_codes{2}, colortables{2}] = read_annotation('/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/fsaverage/label/rh.HCP-MMP1.annot');

hemis = {'lh','rh'};
contrasts = {'SMC', 'WM'};

%% Loop through hemis/modalities and create new annotation files

for hh = 1:length(hemis)
    hemi = hemis{hh};

    % Combine labels
    combined_label = [];
    for ll = 1:length(ROI_names)
        ROI = ROI_names{ll};
        label = readtable([label_dir hemi '.' ROI '_prob_thresh5.label'], 'FileType','text');
        verts = label.Var1;
        combined_label(end+1:end+length(verts)) = verts;
    end
    combined_label_inds = unique(combined_label)+1;

    label_code = label_codes{hh};
    colortable = colortables{hh};

    % Get annot labels overlapped with probmap union
    labels_in_probmap = label_code(combined_label_inds);
    labels_use = unique(labels_in_probmap);

    % Create annotation file with just overlapped labels
    label_code_out = label_code;
    label_code_out(~ismember(label_code_out, labels_use)) = 0;
    colortable_out = colortable;
    colortable_out.numEntries = length(labels_use)+1;
    bad_label_inds = ~ismember(colortable_out.table(:,5),labels_use);
    colortable_out.table(bad_label_inds,:) = [];
    colortable_out.table(end+1,:) = [255 255 255 0 0];
    colortable_out.struct_names = colortable_out.struct_names(~bad_label_inds);
    colortable_out.struct_names{end+1} = 'unlabeled';
    write_annotation([hemi '.GlasserInProb.annot'], vertices, label_code_out, colortable_out);
end

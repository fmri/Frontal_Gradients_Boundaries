%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to create a large frontal label (ROI)
%%% that encompasses the probabilistic areas that respond to sensory stimulation
%%% and WM task for any of 3 modalities
%%%
%%% Tom Possidente - May 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/');
ccc;

%% Initialize Key Variables
hemis = {'lh', 'rh'};

path_base = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/ROIs/';
lobe_annotation_base = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/fsaverage/label/';

%contrasts = {'aA-aP', 'tA-tP', 'vA-vP', 'f-aP', 'f-tP', 'f-vP'};
contrasts = {'vA-vP','aA-aP'};
N_contrasts = length(contrasts);
label_combined = {[],[]};

for hh = 1:2
    hemi = hemis{hh};
    for cc = 1:N_contrasts
        contrast = contrasts{cc};
        path = [path_base hemi '_' contrast '_cortex_probabilistic_thresh5.label'];
        label = readtable(path, 'FileType','text');
        label_combined{hh} = [label_combined{hh}; label];
    end
    label_combined{hh} = unique(label_combined{hh}(:,1:4), 'rows');
    label_combined{hh}(:,5) = table(zeros(height(label_combined{hh}),1)); % last column should be all 0s
    
    %% Load frontal lobe annotation
    lobe_annotation_path = [lobe_annotation_base hemi '.PALS_B12_Lobes.annot'];
    [verts, labels, ctable] = read_annotation(lobe_annotation_path);
    frontal_code = ctable.table(ismember(ctable.struct_names, 'LOBE.FRONTAL'), 5);
    frontal_verts = find(labels==frontal_code)-1; % annots and labels are off by one so subtract 1

    %% Get vertices in combined label that are also in frontal lobe
    frontal_label_mask = ismember(label_combined{hh}{:,1}, frontal_verts);
    label_combined{hh} = label_combined{hh}(frontal_label_mask,:)

    %% Create new label file
    label_fname = ['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/' hemi '.frontal_VisAud_WM_combined.label'];
    label_file = fopen(label_fname,'w');
    fprintf(label_file, ['#!ascii label  , from subject  vox2ras=TkReg\n' num2str(size(label_combined{hh},1)) '\n']);
    writetable(label_combined{hh}, label_fname, 'Delimiter', 'tab', 'WriteMode', 'append', 'FileType', 'text');
    fclose(label_file);
end











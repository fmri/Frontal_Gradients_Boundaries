%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to create frontal sensory and WM labels
%%% for later RAS coordinate analysis. We do this by taking the full cortex
%%% sensory and WM labels and cutting them to just frontal cortex. 
%%%
%%% Tom Possidente - May 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/');
ccc;

%% Initialize Key Variables
hemis = {'lh', 'rh'};

path_base = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/intersections/';
lobe_annotation_base = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/fsaverage/label/';

contrasts = {'f-vP+f-aP', 'vA-vP+aA-aP', 'f-vP+f-aP+f-tP', 'vA-vP+aA-aP+tA-tP'};
N_contrasts = length(contrasts);

for hh = 1:2
    hemi = hemis{hh};

    %% Load frontal lobe annotation
    lobe_annotation_path = [lobe_annotation_base hemi '.PALS_B12_Lobes.annot'];
    [verts, labels, ctable] = read_annotation(lobe_annotation_path);
    frontal_code = ctable.table(ismember(ctable.struct_names, 'LOBE.FRONTAL'), 5);
    frontal_verts = find(labels==frontal_code)-1; % annots and labels are off by one so subtract 1

    for cc = 1:N_contrasts
        contrast = contrasts{cc};
        path = [path_base hemi '_' contrast '_intersection_probabilistic_thresh5.label'];
        label_new = readtable(path, 'FileType','text');

        %% Get vertices in combined label that are also in frontal lobe
        frontal_label_mask = ismember(label_new{:,1}, frontal_verts);
        label_new = label_new(frontal_label_mask,:);

        %% Create new label file
        label_fname = ['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/' hemi '.' contrast '_frontal_probabilistic_thresh5.label'];
        label_file = fopen(label_fname,'w');
        fprintf(label_file, ['#!ascii label  , from subject  vox2ras=TkReg\n' num2str(size(label_new,1)) '\n']);
        writetable(label_new, label_fname, 'Delimiter', 'tab', 'WriteMode', 'append', 'FileType', 'text');
        fclose(label_file);
    end


end








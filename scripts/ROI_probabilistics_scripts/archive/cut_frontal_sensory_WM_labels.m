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

path_base = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/';

%contrasts = {'f-vP+f-aP', 'vA-vP+aA-aP', 'f-vP+f-aP+f-tP', 'vA-vP+aA-aP+tA-tP'};
contrasts = {'f-aA', 'f-vA'};
N_contrasts = length(contrasts);
thresh = '4';
%region = {'lateral', 'medial', 'precentral_sulcus_gyrus'};
region = {'lateral', 'medial'};

for hh = 1:length(hemis)
    hemi = hemis{hh};
    for mm = 1:length(region) % lateral, medial
        %% Load frontal lobe annotation
        frontal_label_path = [path_base hemi '.frontal_' region{mm} '_cortex.label'];
        frontal_label = readtable(frontal_label_path, 'FileType','text');
        frontal_verts = frontal_label.Var1; 

        for cc = 1:N_contrasts
            contrast = contrasts{cc};
            path = [path_base 'intersections/' hemi '_' contrast '_intersection_probabilistic_thresh' thresh '.label'];
            label_new = readtable(path, 'FileType','text');

            %% Get vertices in combined label that are also in frontal lobe
            frontal_label_mask = ismember(label_new{:,1}, frontal_verts);
            label_new = label_new(frontal_label_mask,:);

            %% Create new label file
            label_fname = ['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/' hemi '.' contrast '_frontal_' region{mm} '_probabilistic_thresh' thresh '.label'];
            label_file = fopen(label_fname,'w');
            fprintf(label_file, ['#!ascii label  , from subject  vox2ras=TkReg\n' num2str(size(label_new,1)) '\n']);
            writetable(label_new, label_fname, 'Delimiter', 'tab', 'WriteMode', 'append', 'FileType', 'text');
            fclose(label_file);
        end
    end


end








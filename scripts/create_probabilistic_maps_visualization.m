%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to create a probabilistic map for
%%% how many subjs pass a statistical threshold on particular contrasts at
%%% each vertex
%%% Tom Possidente - July 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/');
ccc;

%% Get Subj Codes

subjCodes{1} = {'MM', 'PP', 'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS', 'AI'};
N{1} = length(subjCodes{1});
subjCodes{2} = {'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS', 'AI'};
N{2} = length(subjCodes{2});


%% Initialize variables
contrasts = {'vA-vP', 'aA-aP'};
contrast_subjs = [1,1];
reverse_contrast = [false false];
intersection = true; % false to loop over contrasts separately, true to make probabilistic of intersection of 2 contrasts

data_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';

p_thresh = 1.3;
N_contrasts = length(contrasts);
hemis = {'lh', 'rh'};
N_hemis = length(hemis);
N_vertices = 163842;


%% Loop through subjs, extract contrast stat data, make probabilistic nii
if intersection
    assert(~any(reverse_contrast), 'reverse contrasts for intersection contrasts not yet implemented')
    assert(contrast_subjs(1)==contrast_subjs(2), 'intersection contrasts with different Ns for each contrast not yet implemented')
    for hh = 1:N_hemis
        hemi = hemis{hh};
        subj_vertices.vol = zeros(1,N_vertices);
        for ss = 1:N{contrast_subjs(1)}
            subjCode = subjCodes{contrast_subjs(1)}{ss};
            stat_path1 = [data_dir subjCode '/localizer/localizer_contrasts_' hemi '/' contrasts{1} '/sig.nii.gz'];
            stat_data1 = MRIread(stat_path1);
            stat_path2 = [data_dir subjCode '/localizer/localizer_contrasts_' hemi '/' contrasts{2} '/sig.nii.gz'];
            stat_data2 = MRIread(stat_path2);

            %%%
            % stat_path3 = [data_dir subjCode '/localizer/localizer_contrasts_' hemi '/vA-aA/sig.nii.gz'];
            % stat_data3 = MRIread(stat_path3);
            %%%

            binarized_stat_data = stat_data1.vol >= p_thresh & stat_data2.vol >= p_thresh;
            %binarized_stat_data = stat_data1.vol >= p_thresh & stat_data2.vol >= p_thresh & abs(stat_data3.vol) < p_thresh;

            subj_vertices.vol = subj_vertices.vol + binarized_stat_data;
        end

        % Save nii file

        MRIwrite(subj_vertices, [hemi '_' contrasts{1} '_' contrasts{2} '_intersection_probabilistic_visualization.nii'])
        disp(['Finishing ' hemi ' ' contrasts]);
    end
else
    for cc = 1:N_contrasts
        contrast = contrasts{cc};
        for hh = 1:N_hemis
            hemi = hemis{hh};
            subj_vertices.vol = zeros(1,N_vertices);
            for ss = 1:N{contrast_subjs(cc)}
                subjCode = subjCodes{contrast_subjs(cc)}{ss};
                stat_path = [data_dir subjCode '/localizer/localizer_contrasts_' hemi '/' contrast '/sig.nii.gz'];
                stat_data = MRIread(stat_path);
                if ismember(contrasts{cc}, {'f-vP', 'f-aP', 'f-tP'}) || reverse_contrast(cc)
                    stat_data.vol = -stat_data.vol; % reverse contrast for interpretability
                end
                binarized_stat_data = stat_data.vol >= p_thresh;
                subj_vertices.vol = subj_vertices.vol + binarized_stat_data;
            end

            % Save nii file
            if reverse_contrast(cc)
                splt_contrast = split(contrast,'-');
                contrast_str = [splt_contrast{2} '-' splt_contrast{1}];
            else
                contrast_str = contrast;
            end

            MRIwrite(subj_vertices, [hemi '_' contrast_str '_probabilistic_visualization.nii'])
            disp(['Finishing ' hemi ' ' contrast_str]);
        end
    end
end

%% Compare maps

ROI_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/';
comp_contrasts = {'vA-aA', 'aA-vA'};
hemis = {'lh', 'rh'};

for hh = 1:2
    hemi = hemis{hh};
    map1 = MRIread([ROI_dir hemi '_' comp_contrasts{1} '_probabilistic_visualization.nii']);
    map1_overthresh = map1.vol>=4;
    map2 = MRIread([ROI_dir hemi '_' comp_contrasts{2} '_probabilistic_visualization.nii']);
    map2_overthresh = map2.vol>=4;
    comp_map.vol = map1.vol-map2.vol;
    comp_map.vol(~map2_overthresh & ~map1_overthresh) = 0;
    MRIwrite(comp_map, [ROI_dir hemi '_' comp_contrasts{1} '_minus_' comp_contrasts{2} '_probabilistic_visualization.nii'])
end


%% Make label files from comparison maps
ROI_dir_old = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/ROIs/';
ROI_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/';
comp_contrasts = {'vA-aA', 'aA-vA'};
hemis = {'lh', 'rh'};

lh_cortex_label = readtable([ROI_dir_old 'lh_inflated_wholecortex.label'], 'FileType','text');
rh_cortex_label = readtable([ROI_dir_old 'rh_inflated_wholecortex.label'], 'FileType','text');
lhrh_cortex_label = {lh_cortex_label, rh_cortex_label};

for hh = 1:2
    hemi = hemis{hh};
    map1 = MRIread([ROI_dir hemi '_' comp_contrasts{1} '_probabilistic_visualization.nii']);
    map1_overthresh = map1.vol>=4;
    map2 = MRIread([ROI_dir hemi '_' comp_contrasts{2} '_probabilistic_visualization.nii']);
    map2_overthresh = map2.vol>=4;
    comp_map.vol = map1.vol-map2.vol;
    comp_map.vol(~map2_overthresh & ~map1_overthresh) = 0;

    cortex_label = lhrh_cortex_label{hh};
    binarized_probabilistic_vis = comp_map.vol >= 4;
    binarized_probabilistic_aud = comp_map.vol <= -4;
    label_inds_vis = find(binarized_probabilistic_vis) - 1;
    label_inds_aud = find(binarized_probabilistic_aud) - 1;
    label_rows_vis = cortex_label(ismember(cortex_label.Var1, label_inds_vis),:);
    label_rows_aud = cortex_label(ismember(cortex_label.Var1, label_inds_aud),:);

    
    % Make label file 
    label_fname = [ROI_dir hemi '_' contrast '_cortex_probabilistic_thresh' num2str(tt) '.label']; % needs to be changed, left off here
    %label = unique_rows(counts>=tt,:);
    label_file = fopen(label_fname,'w');
    fprintf(label_file, ['#!ascii label  , from subject  vox2ras=TkReg\n' num2str(size(label,1)) '\n']);
    writematrix(label, label_fname, 'Delimiter', 'tab', 'WriteMode', 'append', 'FileType', 'text');
    fclose(label_file);
end





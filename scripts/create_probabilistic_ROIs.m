%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to create probabilistic ROIs from the
%%% passive-fixation (sensory drive) and active-passive (working memory) contrasts
%%% of subjs localizer data by finding areas of overlap
%%% Tom Possidente - August 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/');
ccc;

%% Get Subj Codes

experiment_name = 'spacetime';

ROI_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/ROIs/';
ROI_dir_save = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/intersections/';

subjCodes{1} = {'MM', 'PP', 'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS', 'AI'};
N{1} = length(subjCodes{1});
subjCodes{2} = {'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS', 'AI'};
N{2} = length(subjCodes{2});


%% Initialize variables
data_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';
label_or_nii = 'nii'; % 'label' to create .label files, 'nii' to create nifti files with # of subjs as data
t_thresh = 2;
% contrasts = {'f-vP', 'f-aP', 'f-tP', 'vA-vP', 'aA-aP', 'tA-tP'};
% contrast_subjs = [2,2,2,1,1,1];
contrasts = {'vA-aA'};
contrast_subjs = [1,1];
N_contrasts = length(contrasts);
hemis = {'lh', 'rh'};
N_hemis = length(hemis);
N_vertices = 163842;
prob_ROI_label = cell(N_contrasts, N_hemis);

lh_cortex_label = readtable([ROI_dir 'lh_inflated_wholecortex.label'], 'FileType','text');
rh_cortex_label = readtable([ROI_dir 'rh_inflated_wholecortex.label'], 'FileType','text');
lhrh_cortex_label = {lh_cortex_label, rh_cortex_label};

%% Loop through subjs and extract contrast tstat data

for cc = 1:N_contrasts
    contrast = contrasts{cc};
    for hh = 1:N_hemis
        subj_vertices = [];
        for ss = 1:N{contrast_subjs(cc)}
            subjCode = subjCodes{contrast_subjs(cc)}{ss};

            cortex_label = lhrh_cortex_label{hh};
            hemi = hemis{hh};
            tstat_path = [data_dir subjCode '/localizer/localizer_contrasts_' hemi '/' contrast '/t.nii.gz'];
            tstat_data = MRIread(tstat_path);
            if ismember(contrasts{cc}, {'f-vP', 'f-aP', 'f-tP'})
                tstat_data.vol = -tstat_data.vol; % reverse contrast for interpretability
            end
            binarized_tstat_data = tstat_data.vol >= t_thresh;
            label_inds = find(binarized_tstat_data) - 1;
            label_rows = cortex_label(ismember(cortex_label.Var1, label_inds),:);
            nrows = size(subj_vertices, 1);
            subj_vertices(nrows+1:nrows+size(label_rows,1),:) = table2array(label_rows);
        end
        [unique_rows,~,ind] = unique(subj_vertices,'rows');
        counts = histc(ind,unique(ind)); % count occurences of each vertex
        unique_rows(:,5) = counts;
        prob_ROI_label{cc,hh} = unique_rows;

        % Save label files with different thresholds
        for tt = 3:max(counts)-1
            label_fname = [ROI_dir hemi '_' contrast '_cortex_probabilistic_thresh' num2str(tt) '.label'];
            label = unique_rows(counts>=tt,:);
            label_file = fopen(label_fname,'w');
            fprintf(label_file, ['#!ascii label  , from subject  vox2ras=TkReg\n' num2str(size(label,1)) '\n']);
            writematrix(label, label_fname, 'Delimiter', 'tab', 'WriteMode', 'append', 'FileType', 'text');
            fclose(label_file);
        end
        disp(['Finishing ' hemi ' ' contrast]);
    end
end

%% Now create probabilistic intersection labels

contrasts = {'vA-vP', 'aA-aP'}; % Contrasts to get probabilistic intersections for
contrast_subjs = 1;
N_contrasts = length(contrasts);
prob_ROI_label = cell(N_hemis,1);

for hh = 1:N_hemis
    hemi = hemis{hh};
    cortex_label = lhrh_cortex_label{hh};
    all_subj_vertices = [];
    for ss = 1:N{contrast_subjs}
        subjCode = subjCodes{contrast_subjs}{ss};
        subj_vertices = [];
        for cc = 1:N_contrasts
            contrast = contrasts{cc};
            if cc==1
                contrast_str = contrast;
            else
                contrast_str = [contrast_str '+' contrast];
            end
            tstat_path = [data_dir subjCode '/localizer/localizer_contrasts_' hemi '/' contrast '/t.nii.gz'];
            tstat_data = MRIread(tstat_path);
            if ismember(contrasts{cc}, {'f-vP', 'f-aP', 'f-tP'})
                tstat_data.vol = -tstat_data.vol; % reverse contrast for interpretability
            end
            binarized_tstat_data = tstat_data.vol >= t_thresh;
            label_inds = find(binarized_tstat_data) - 1;
            label_rows = cortex_label(ismember(cortex_label.Var1, label_inds),:);
            nrows = size(subj_vertices, 1);
            subj_vertices(nrows+1:nrows+size(label_rows,1),:) = table2array(label_rows);
        end
        [unique_rows,~,ind] = unique(subj_vertices,'rows');
        counts = histc(ind,unique(ind)); % count occurences of each vertex
        unique_rows = unique_rows(counts==N_contrasts,:);
        all_subj_vertices = [all_subj_vertices; unique_rows];
    end
    [unique_rows,~,ind] = unique(all_subj_vertices,'rows');
    counts = histc(ind,unique(ind)); % count occurences of each vertex
    unique_rows(:,5) = counts;
    prob_ROI_label{hh} = unique_rows;

    % Save label files with different thresholds
    for tt = 3:max(counts)-1
        label_fname = [ROI_dir_save hemi '_' contrast_str '_intersection_probabilistic_thresh' num2str(tt) '.label'];
        label = unique_rows(counts>=tt,:);
        label_file = fopen(label_fname,'w');
        fprintf(label_file, ['#!ascii label  , from subject  vox2ras=TkReg\n' num2str(size(label,1)) '\n']);
        writematrix(label, label_fname, 'Delimiter', 'tab', 'WriteMode', 'append', 'FileType', 'text');
        fclose(label_file);
    end
    disp(['Finished ' hemi ' ' contrast_str]);
end

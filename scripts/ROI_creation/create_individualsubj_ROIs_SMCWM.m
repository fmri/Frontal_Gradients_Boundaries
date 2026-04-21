%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to create subject-specific ROIs using
%%% probabilistic ROI search spaces and individual contrast data
%%%
%%% Tom Possidente - January 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables
ROIs = {'aINS', 'aIPS', 'aMFG', 'ant_temporal', 'cIPS', 'DO', 'inf_lat_frontal', 'LOT', 'midIFS', 'midINS', ...
    'ms_post_STSG', 'MT', 'parietal_opercular', 'pIPS', 'post_col_sulc', 'post_temporal', 'preSMA', 'sup_lat_frontal',...
    'VO', 'VOT'};
N_ROIs = length(ROIs);

hemis = {'lh', 'rh'};

subjCodes = {'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS', 'AI'};
N = length(subjCodes);

contrasts = {'vAvP-f', 'aAaP-f'};
N_contrasts = length(contrasts);
contrast_names = {'visual', 'auditory', 'supramodal'};

N_vertices = 163842;
vertex_inds = 1:N_vertices;

pval_thresh = -log10(0.01); % pvals in freesurfer files are expressed as -log10(p)
ROI_vertices_thresh = 0.1; % percent of vertices in search space that must exist in final ROI
ROI_subj_thresh = 0.5; % precent of subjs that must have the ROI for it to be used

probROI_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/probabilistic_allROIs/';
data_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';
ROI_outdir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/subj_specific_01/';

% Load cortex labels for later use
cortex_label_lh = readtable('/projectnb/somerslab/tom/fs_helpful_files/lh_cortex.label', 'FileType','text');
cortex_label_rh = readtable('/projectnb/somerslab/tom/fs_helpful_files/rh_cortex.label', 'FileType','text');
cortex_labels = {cortex_label_lh, cortex_label_rh};

%% Preload probabilistic ROI masks
probabilistic_ROIs = nan(N_vertices, N_ROIs, 2);
for hh = 1:2
    hemi = hemis{hh};
    for rr = 1:N_ROIs
        ROI = ROIs{rr};

        % Preload probabilistic label
        label = readtable([probROI_dir hemi '.' ROI '_prob_thresh5.label'], 'FileType','text');
        vert_mask = ismember(vertex_inds, label{:,1}+1);
        probabilistic_ROIs(:,rr,hh) = vert_mask;
        disp([hemi ' ' ROI ': ' num2str(sum(vert_mask)) ' vertices'])
    end
end

%% Loop through subjs/ROIs/contrasts to create ROIs
ROI_size = nan(N, 2, N_ROIs, N_contrasts+1); % subjs x hemis x ROIs x contrast+1
ROI_good = nan(N, 2, N_ROIs, N_contrasts+1);

for ss = 1:N
    subjCode = subjCodes{ss};
    for hh = 1:2
        hemi = hemis{hh};
        for rr = 1:N_ROIs
            ROI_name = ROIs{rr};
            prob_ROI_mask = probabilistic_ROIs(:,rr,hh);
            final_ROImasks = nan(N_vertices, N_contrasts);
            for cc = 1:N_contrasts
                contrast = contrasts{cc};
                data = MRIread([data_dir subjCode '/localizer/localizer_contrasts_' hemi '/' contrast '/sig.nii.gz']);
                data_thresh_mask = data.vol >= pval_thresh;
                final_ROImask = prob_ROI_mask & data_thresh_mask';
                final_ROImasks(:,cc) = final_ROImask;
                ROI_size(ss,hh,rr,cc) = sum(final_ROImask);
                ROI_good(ss,hh,rr,cc) = sum(final_ROImask) >= (sum(prob_ROI_mask)*ROI_vertices_thresh); % are the number of ROIs more than 10% of the search space?

                % Create ROI
                if ROI_good(ss,hh,rr,cc)
                    cortex_label_inds = vertex_inds(final_ROImask) - 1;
                    cortex_label_mask = ismember(cortex_labels{hh}{:,1}, cortex_label_inds);
                    assert(length(cortex_label_inds)==sum(cortex_label_mask), 'number of vertices in ROI mask does not match number of vertices in label');
                    label = table2array(cortex_labels{hh}(cortex_label_mask,:));
                    label_fname = [ROI_outdir hemi '.' ROI_name '_' contrast_names{cc} '_' subjCode '.label'];
                    if ~isfile(label_fname)
                        label_file = fopen(label_fname,'w');
                        fprintf(label_file, ['#!ascii label  , from subject  vox2ras=TkReg\n' num2str(size(label,1)) '\n']);
                        writematrix(label, label_fname, 'Delimiter', 'tab', 'WriteMode', 'append', 'FileType', 'text');
                        fclose(label_file);
                    end
                end
            end

            % Create intersection of 2 contrast ROIs
            final_ROImask_intersection = final_ROImasks(:,1) & final_ROImasks(:,2);
            ROI_size(ss,hh,rr,N_contrasts+1) = sum(final_ROImask_intersection);
            ROI_good(ss,hh,rr,N_contrasts+1) = sum(final_ROImask_intersection) >= (sum(prob_ROI_mask)*0.1);

            % Create supramodal ROI
            if ROI_good(ss,hh,rr,3)
                cortex_label_inds = vertex_inds(final_ROImask_intersection) - 1;
                cortex_label_mask = ismember(cortex_labels{hh}{:,1}, cortex_label_inds);
                assert(length(cortex_label_inds)==sum(cortex_label_mask), 'number of vertices in ROI mask does not match number of vertices in label');
                label = table2array(cortex_labels{hh}(cortex_label_mask,:));
                label_fname = [ROI_outdir hemi '.' ROI_name '_' contrast_names{3} '_' subjCode '.label'];
                if ~isfile(label_fname)
                    label_file = fopen(label_fname,'w');
                    fprintf(label_file, ['#!ascii label  , from subject  vox2ras=TkReg\n' num2str(size(label,1)) '\n']);
                    writematrix(label, label_fname, 'Delimiter', 'tab', 'WriteMode', 'append', 'FileType', 'text');
                    fclose(label_file);
                end
            end
        end
    end
    disp(['finished subj ' subjCode]);
end

%% Check %s of ROIs
perc_ROIs = nan(N_ROIs, 2, N_contrasts+1);
for rr = 1:N_ROIs
    for hh = 1:2
        for cc = 1:N_contrasts+1
            perc_ROIs(rr,hh,cc) = mean(ROI_good(:,hh,rr,cc));
        end
    end
end

%% Plot ROIs Ts
contrast_names = {'visual', 'auditory', 'supramodal'};
for cc = 1:N_contrasts+1
    figure;
    bar(1:N_ROIs, [perc_ROIs(:,1,cc), perc_ROIs(:,2,cc)]);
    xticks(1:N_ROIs);
    xticklabels(strrep(ROIs, '_', ' '));
    legend({'lh', 'rh'});
    ylabel('Percent Good ROIs');
    xlabel('ROI');
    title([contrast_names{cc} ' | ' num2str(round(mean(perc_ROIs(:,:,cc), 'all'),3) )])
    yline(.50, '--r');
end


%% Reject ROIs with less than 50% subjs
contrast_names = {'visual', 'auditory', 'supramodal'};
for cc = 1:N_contrasts+1
    good_ROIs = sum([perc_ROIs(:,1,cc), perc_ROIs(:,2,cc)]>ROI_subj_thresh,2)==2; % both hemis have over 50%

    % Delete ROI labels for ROIs without enough subjs
    for rr = 1:N_ROIs
        if ~good_ROIs(rr)
            disp(['rm ' ROI_outdir 'lh.' ROIs{rr} '_' contrast_names{cc} '_*'])
            disp(['rm ' ROI_outdir 'rh.' ROIs{rr} '_' contrast_names{cc} '_*'])

            unix(['rm ' ROI_outdir 'lh.' ROIs{rr} '_' contrast_names{cc} '_*'])
            unix(['rm ' ROI_outdir 'rh.' ROIs{rr} '_' contrast_names{cc} '_*'])
        end
    end

    new_perc_data_ROI = [perc_ROIs(good_ROIs,1,cc), perc_ROIs(good_ROIs,2,cc)];
    disp([contrast_names{cc} ' ROIs removed: ' ROIs(~good_ROIs)]);

    figure;
    bar(1:sum(good_ROIs), new_perc_data_ROI);
    xticks(1:sum(good_ROIs));
    xticklabels(strrep(ROIs(good_ROIs), '_', ' '));
    legend({'lh', 'rh'});
    ylabel('Percent Subjs with ROI');
    xlabel('ROI');
    title([contrast_names{cc} ' | ' num2str(round(mean(new_perc_data_ROI, 'all'),3) )])
    yline(.50, '--r');
end
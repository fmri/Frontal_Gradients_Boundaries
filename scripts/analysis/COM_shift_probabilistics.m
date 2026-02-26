%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to calculate the center of mass of
%%% permuted probabilistic ROIs for WM and SMC contrasts, then compare those flatmap
%%% XY locations (direction and distance)
%%%
%%% Tom Possidente - January 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables
perm_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/permutation_probabilistic_maps_WMSMC/';
ROI_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/probabilistic_allROIs/';
unpermuted_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/SMC_WM_nonpermuted_probabilistics/';

hemis = {'lh', 'rh'};
contrasts = {'WM', 'SMC'};
permutation_names = {'group1', 'group2'};
keyword = 'visual';

% Establish which ROIs are included in the analysis (based on mean PSC ratio threshold and number of subjs with the ROI
switch keyword
    case {'visual', 'vis'}
        ROI_names = {'aMFG', 'midIFS', 'aINS', 'preSMA', 'inf_lat_frontal', 'sup_lat_frontal', ...
                     'aIPS', 'VOT', 'cIPS', 'LOT', 'pIPS', 'VO', 'DO'};
    case {'auditory', 'aud'}
        ROI_names = {'aMFG', 'midIFS', 'aINS', 'preSMA', 'inf_lat_frontal', 'sup_lat_frontal', ...
                     'parietal_opercular', 'aIPS', 'ms_post_STSG', 'cIPS'};
    case {'supramodal', 'supra'}
        ROI_names = {'aINS', 'preSMA', 'inf_lat_frontal', 'sup_lat_frontal', 'aIPS', 'cIPS', 'midIFS'};
end


N_iterations = 10000;
N_hemis = length(hemis);
N_ROIs = length(ROI_names);
N_contrasts = length(contrasts);

N_thresh = 5; % 5 subj threshold for probabilistic map

COMs = nan(N_ROIs, N_contrasts, N_hemis, N_iterations, 2); % storage for center of masses XY coordinates

plotting_diagnostics = false;

%% Loop through ROI probabilistic data and calculate center of mass for each
tic;

for cc = 1:N_contrasts
    group = permutation_names{cc};
    for rr = 1:N_ROIs
        ROI_name = ROI_names{rr};
        for hh = 1:N_hemis
            hemi = hemis{hh};

            % Load patch file
            patch_path = [ROI_dir hemi '.' ROI_name '_prob_thresh5_flat.patch'];
            patch = read_patch(patch_path);

            % Load label file
            label_path = [ROI_dir hemi '.' ROI_name '_prob_thresh5.label'];
            label = readtable(label_path, 'FileType', 'text');

            % Get #subj in vertices in both patch and label (there may be some negligable differences due to imperfect cutting
            inds_label = label{:,1};
            patch_in_label = ismember(inds_label, patch.ind);
            perc_labelinpatch = sum(patch_in_label)/length(inds_label);
            if perc_labelinpatch<0.99
                warning([hemi ' ' ROI_name ' ' group ' percent of label in patch was ' num2str(perc_labelinpatch)])
            end
            inds = inds_label(patch_in_label);

            for ii = 1:N_iterations
                % Load probabilistic data for full hemisphere
                prob_data_path = [perm_dir hemi '_iter' num2str(ii) '_' group '_' keyword '.nii'];
                prob_data = MRIread(prob_data_path);

                % Restrict to prob map >= N_thresh
                ROI_inds_good = find(prob_data.vol >= N_thresh);
                inds_iter = inds(ismember(inds, ROI_inds_good-1)); % -1 to convert from nii indexing to label indexing
                inds_iter = sort(inds_iter); % keep order consistent
                ROI_data = prob_data.vol(inds_iter+1); % +1 to convert from label indexing to nii indexing

                % Calculate center of mass using XY coords and #subjs per vertex
                ind_mask = ismember(patch.ind, inds_iter);
                x = patch.x(ind_mask);
                y = patch.y(ind_mask);
                patch_inds = patch.ind(ind_mask);
                COMs(rr,cc,hh,ii,:) = ROI_data * [x; y]' / sum(ROI_data);

                if plotting_diagnostics
                    % if length(findobj('type','figure'))>50
                    %     disp('more than 50 plots already open, skipping plotting');
                    %     continue;
                    % end
                    % 
                    % figure;
                    % scatter(x,y,[],ROI_data, 'filled'); hold on;
                    % scatter(COMs(rr,cc,hh,ii,1), COMs(rr,cc,hh,ii,2), 100, 'r', 'filled');
                    % %scatter(x(closest_ind), y(closest_ind), 100, 'o', 'filled');
                    % title([hemi ' ' ROI_name ' ' group ' iter' num2str(ii)])
                end

            end
        end
    end
end

dimensions = {'ROI', 'group', 'hemisphere', 'iteration', 'XYcoordinates'};
save(['permutation_COMs_' keyword '_' num2str(N_iterations) 'iters.mat'], 'COMs', 'dimensions');
toc;

%% Calculate COM of original WM and SMC grouped data
original_COMs = nan(N_ROIs, N_contrasts, N_hemis, 2);
for hh = 1:N_hemis
    hemi = hemis{hh};
    for cc = 1:N_contrasts
        group = permutation_names{cc};
        contrast = contrasts{cc};

        % Load probabilistic data for full hemisphere
        prob_data_path = [perm_dir hemi '_original_' contrast '_' keyword '.nii'];
        prob_data = MRIread(prob_data_path);

        ROI_COM_label_data = nan(N_ROIs,5);

        for rr = 1:N_ROIs
            ROI_name = ROI_names{rr};

            % Load patch file
            patch_path = [ROI_dir hemi '.' ROI_name '_prob_thresh5_flat.patch'];
            patch = read_patch(patch_path);

            % Load label file
            label_path = [ROI_dir hemi '.' ROI_name '_prob_thresh5.label'];
            label = readtable(label_path, 'FileType', 'text');

            % Get #subj in vertices in both patch and label (there will be some negligable differences due to imperfect cutting
            inds_label = label{:,1};
            patch_in_label = ismember(inds_label, patch.ind);
            perc_labelinpatch = sum(patch_in_label)/length(inds_label);
            if perc_labelinpatch<0.99
                warning([hemi ' ' ROI_name ' ' group ': percent of label in patch was ' num2str(perc_labelinpatch)])
            end
            inds = inds_label(patch_in_label);

            % Restrict to prob map >= N_thresh
            ROI_inds_good = find(prob_data.vol >= N_thresh);
            inds = inds(ismember(inds, ROI_inds_good-1)); % -1 to convert from nii indexing to label indexing
            inds = sort(inds); % keep order consistent
            ROI_data = prob_data.vol(inds+1); % +1 to convert from label indexing to nii indexing

            % Calculate center of mass using XY coords and #subjs per vertex
            ind_mask = ismember(patch.ind, inds);
            x = patch.x(ind_mask);
            y = patch.y(ind_mask);
            patch_inds = patch.ind(ind_mask);
            original_COMs(rr,cc,hh,:) = ROI_data * [x; y]' / sum(ROI_data);

            % Find and record vertex index closest to COM point
            if ~any(isnan(original_COMs(rr,cc,hh,:)))
                [vert_dist, closest_ind] = min( sqrt( sum( ([x;y] - squeeze(original_COMs(rr,cc,hh,:))).^2 ) ) ); 
                if vert_dist>0.999
                    keyboard;
                end
                closest_label_ind = patch_inds(closest_ind);
                ROI_COM_label_data(rr,:) = label{label.Var1==closest_label_ind,:};
            end

            if plotting_diagnostics
                if length(findobj('type','figure'))>100
                    disp('more than 100 plots already open, skipping plotting');
                    continue;
                end

                figure;
                scatter(x,y,[],ROI_data, 'filled'); hold on;
                scatter(COMs(rr,cc,hh,1), COMs(rr,cc,hh,ii,2), 100, 'r', 'filled');
                %scatter(x(closest_ind), y(closest_ind), 100, 'o', 'filled');
                title([hemi ' ' ROI_name ' ' contrast])
            end

        end

        % Create label with COM points for all ROIs
        label_fname = [unpermuted_dir hemi '.' keyword '_' contrast '_COMs.label'];
        label_file = fopen(label_fname,'w');
        fprintf(label_file, ['#!ascii label  , from subject  vox2ras=TkReg\n' num2str(size(ROI_COM_label_data,1)) '\n']);
        writetable(array2table(ROI_COM_label_data), label_fname, 'Delimiter', 'tab', 'WriteMode', 'append', 'FileType', 'text');
        fclose(label_file);
    end
end

%% Compare direction and distance between COMs
tic;
coord_diffs = table();

counter = 0;
for ii = 1:N_iterations
    for hh = 1:N_hemis
        hemi = hemis{hh};
        for rr = 1:N_ROIs
            ROI_name = ROI_names{rr};
            counter = counter + 1;
            COM_g1 = squeeze(COMs(rr,1,hh,ii,:));
            COM_g2 = squeeze(COMs(rr,2,hh,ii,:));

            if plotting_diagnostics & strcmp(ROI_name, 'iPCS')
                patch_path = [ROI_dir hemi '.' ROI_name '_prob_thresh5_flat.patch'];
                patch = read_patch(patch_path);
                figure;
                scatter(patch.x, patch.y); hold on;
                scatter(COM_g1(1), COM_g1(2), 25, 'filled');
                scatter(COM_g2(1), COM_g2(2), 25, 'filled');
                legend({'coordinates', 'WM peak', 'SMC peak'});
                title([hemis{hh} ' ' ROI_name ' iter' num2str(ii) ' dist=' num2str(round(pdist([COM_g1'; COM_g2']), 2))]);
            end
            
            add_table = table(hemis(hh), {ROI_name}, ii, (COM_g1-COM_g2)', pdist([COM_g1'; COM_g2']));
            coord_diffs = [coord_diffs; add_table];
        end
    end
end
coord_diffs.Properties.VariableNames = {'hemisphere', 'ROI', 'iteration', 'XY_difference', 'distance'};
toc;

save(['permutation_COM_distances_' keyword '_' num2str(N_iterations) 'iters.mat'], 'coord_diffs');

%% Calculate distance for original WM and SMC grouped data
orig_coord_diffs = table();

counter = 0;
for hh = 1:N_hemis
    hemi = hemis{hh};
    for rr = 1:N_ROIs
        ROI_name = ROI_names{rr};
        counter = counter + 1;
        COM_WM = squeeze(original_COMs(rr,1,hh,:));
        COM_SMC = squeeze(original_COMs(rr,2,hh,:));

        patch_path = [ROI_dir hemi '.' ROI_name '_prob_thresh5_flat.patch'];
        patch = read_patch(patch_path);

        if plotting_diagnostics
            figure;
            scatter(patch.x, patch.y); hold on;
            scatter(COM_WM(1), COM_WM(2), 25, 'filled');
            scatter(COM_SMC(1), COM_SMC(2), 25, 'filled');
            legend({'coordinates', 'WM peak', 'SMC peak'});
            title([hemis{hh} ' ' ROI_name]);
        end

        orig_coord_diffs{counter,:} = {hemis{hh}, ROI_name, squeeze(COM_WM-COM_SMC)', pdist([COM_WM'; COM_SMC'])};
    end
end
orig_coord_diffs.Properties.VariableNames = {'hemisphere', 'ROI', 'XY_difference', 'distance'};

%% Plot null distriubtion of distances
results_table = table();
for hh = 1:N_hemis
    hemi = hemis{hh};
    for rr = 1:N_ROIs
        ROI_name = ROI_names{rr};
        perm_data = coord_diffs.distance(ismember(coord_diffs.hemisphere, hemi) & ismember(coord_diffs.ROI, ROI_name));
        orig_data = orig_coord_diffs.distance{ismember(orig_coord_diffs.hemisphere, hemi) & ismember(orig_coord_diffs.ROI, ROI_name)};

        figure;
        histogram(perm_data, 'NumBins',100); hold on;
        xline(orig_data, '--r');
        title([hemi ' ' ROI_name ' ' keyword ' permutation p-value=' num2str(sum(orig_data<perm_data)/(N_iterations+1))])
        
        pval = sum(orig_data<perm_data)/(N_iterations+1);
        disp([hemi ' ' ROI_name ' permutation p-value=' num2str(pval)]);
        results_table = [results_table; {keyword, hemi, ROI_name, orig_data, pval}];
    end
end
results_table.Properties.VariableNames = {'Modality', 'Hemisphere', 'ROI', 'COM_dist', 'perm_pval'};
save([keyword '_permtest_results.mat'], 'results_table')

%% FDR testing

load('/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/permutation_analysis_results/auditory_permtest_results.mat','results_table');
results_auditory = results_table;
load('/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/permutation_analysis_results/visual_permtest_results.mat','results_table');
results_visual = results_table;
load('/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/permutation_analysis_results/supramodal_permtest_results.mat','results_table');
results_supramodal = results_table;

results_table_all = [results_auditory; results_visual; results_supramodal];
[c_p,c_a,h] = fdr_BH(results_table_all.perm_pval, .05);
results_table_all.FDR_rejectnull = h;
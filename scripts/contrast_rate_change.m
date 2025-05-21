%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to get a measure of the rate of change of
%%% different contrasts across the cortex in order to see if there are
%%% abrupt or gradual changes in funcitonal properties
%%% Tom Possidente - May 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/');
ccc;

%% Initialize Key Variables
%reject_subjs = {'AH', 'SL', 'RR'};
subjCodes = {'MM', 'PP', 'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS', 'AI'};
N = length(subjCodes);
hemis = {'lh', 'rh'};

contrast_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';
contrast_name = 'A-P';
gradient_direction = 2; % 1 = Right->Left for lh, opposite for rh, 2 = Posterior->Anterior, 3 = Inferior->Superior 
gradient_dir_strs = {'RL', 'PA', 'IS'};
spotlight = 20; % number of (closest neighbor) vertices used to calculate activation gradient
dist_thresh = 5; % how close (in RAS coords) all neighbors must be to use vertex in gradient calc
gradient_means_SDs = nan(N, 2, 2); 
gradient_absmeans_SDs = nan(N, 2, 2);

%% Load frontal labels 
lh_frontal_label = readtable('/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/lh.frontal_sensory_WM_3modality_combined.label', 'FileType','text');
rh_frontal_label = readtable('/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/rh.frontal_sensory_WM_3modality_combined.label', 'FileType','text');
frontal_labels = {lh_frontal_label, rh_frontal_label};

%% Precalculate distances between vertices and neighbors in label (saves lots of computing time)
% vert_neighbors = {{},{}};
% for hh = 1:2
%     for vv = 1:height(frontal_labels{hh})
%         vert = frontal_labels{hh}{vv,:};
%         distances = sum( (vert(2:4) - frontal_labels{hh}{:,2:4}).^2 , 2); % distances between all vertices and vert
%         [sorted_dists, dist_inds] = sort(distances);
%         if any(sorted_dists(1:spotlight)>dist_thresh)
%             vert_neighbors{hh}{vv} = nan;
%         else
%             vert_neighbors{hh}{vv} = dist_inds(1:spotlight);
%         end
%     end
% end
% save("frontal_label_vert_50neighbors_10thresh.mat", 'vert_neighbors');

load(['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/frontal_label_vert_neighbors' num2str(spotlight) '_' num2str(dist_thresh) 'thresh.mat'], 'vert_neighbors');

%% Loop through subjs and compute cortical gradient

for ss = 1:N
    subjCode = subjCodes{ss};
    for hh = 1:2
        count = 0;
        hemi = hemis{hh};

        %% Load t stats for contrast
        path = [contrast_dir subjCode '/localizer/localizer_contrasts_' hemi '/' contrast_name '/t.nii.gz'];
        data = MRIread(path);

        %% Get tstats front frontal label only
        label_data = data.vol(frontal_labels{hh}{:,1}+1); % label file vertex indices are off by one, so add 1
        frontal_labels{hh}{:,5} = label_data';

        %% Loop through vertices, extract neighbors, extract gradient
        for vv = 1:height(frontal_labels{hh})       
            if isnan(vert_neighbors{hh}{vv})
                frontal_labels{hh}{vv,6} = nan;
                count = count + 1;
            else
                neighbors = frontal_labels{hh}(vert_neighbors{hh}{vv},:);
                fit = polyfit(neighbors{:,gradient_direction+1}, neighbors{:,5}, 1); % might be able to find a faster function to calc slope here if speedup is needed
                gradient = fit(1);
                frontal_labels{hh}{vv,6} = gradient;
                % neighbors = frontal_labels{hh}(vert_neighbors{hh}{vv},:);
                % fit_A = polyfit(neighbors{:,3}, neighbors{:,5}, 1); 
                % fit_S = polyfit(neighbors{:,4}, neighbors{:,5}, 1); 
                % angle = atan2d(fit_A(1), fit_S(1));
                % mag = sqrt(fit_A(1)^2 + fit_S(1)^2);
                % if mag >= 0.05
                %     frontal_labels{hh}{vv,6} = angle;
                % else
                %     frontal_labels{hh}{vv,6} = nan;
                %     count = count + 1;
                % end
            end
        end
        disp(['Percentage vertices without enough close neighbors: ' num2str(count/height(frontal_labels{hh}))])
        %disp(['Percentage vertices without enough close neighbors or weak gradient: ' num2str(count/height(frontal_labels{hh}))])

        % Get means/SDs of gradients
        gradient_means_SDs(ss,hh,1) = mean(frontal_labels{hh}{:,6}, 'omitnan');
        gradient_means_SDs(ss,hh,2) = std(frontal_labels{hh}{:,6}, 'omitnan');
        gradient_absmeans_SDs(ss,hh,1) = mean(abs(frontal_labels{hh}{:,6}), 'omitnan');
        gradient_absmeans_SDs(ss,hh,2) = std(abs(frontal_labels{hh}{:,6}), 'omitnan');

        %% Write gradient nii
        outfile = ['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/gradient_niftis/' hemi '_' subjCode '_' contrast_name '_' gradient_dir_strs{gradient_direction} '_' num2str(spotlight) '_gradient.nii.gz'];
        %outfile = ['/projectnb/somerslab/tom/projects/sensory_WM_topology/data/gradient_niftis/' hemi '_' subjCode '_' contrast_name '_' num2str(spotlight) '_gradient_angles.nii.gz'];
        data.vol = nan(size(data.vol));
        data.vol(frontal_labels{hh}{:,1}+1) = frontal_labels{hh}{:,6};
        MRIwrite(data, outfile)
    end
end

save(['gradient_meanSDs_' contrast_name '_' gradient_dir_strs{gradient_direction} '.mat'], "gradient_means_SDs", "gradient_absmeans_SDs", "subjCodes");





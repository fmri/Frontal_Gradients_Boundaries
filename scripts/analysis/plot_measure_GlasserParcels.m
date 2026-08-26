%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to visualize flattened Glasser parcels
%%% and calculate their distance across
%%%
%%% Tom Possidente - Septemeber 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables 
flatpatch_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/fsaverage/surf/';
parcels = {'IFJp', 'i6-8', 'AVI'};
hemis = {'L', 'R'};

%% plot
for hh = 1:length(hemis)
    hemi = hemis{hh};
    for pp = 1:length(parcels)
        parcel = parcels{pp};
        flatpatch = read_patch([flatpatch_dir hemi '_' parcel '_ROI_flatpatch']);
        x = flatpatch.x - mean(flatpatch.x);
        y = flatpatch.y - mean(flatpatch.y);
        [coef, score, latent] = pca([x;y]', 'NumComponents',2);
        figure;
        scatter(x, y);
        hold on;
        scatter(score(:,1), score(:,2), 'filled')
        PC1_dist = max(score(:,1)) - min(score(:,1))
        PC2_dist = max(score(:,2)) - min(score(:,2))
    end
end
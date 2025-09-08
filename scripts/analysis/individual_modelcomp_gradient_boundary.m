%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to compare a gradient-based and
%%% boundary-based model of change in function (contrast PSC) across an ROI in
%%% individual subjects
%%%
%%% Tom Possidente - Septemeber 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/');
ccc;

%% Initialize Key Variables
contrasts = {'aPvP-f', 'vAaA-vPaP'};
latmed = 'medial'; % lateral or medial

hemis = {'lh', 'rh'};
fs_num = 163842;

experiment_name = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode(~ismember(subjDf_cut.subjCode, {'AH', 'SL', 'RR', 'PP', 'MM'}));
N_subjs = length(subjCodes);

subjdir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';
fs_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/fsaverage/surf/';

axis_method = 'average'; % regression or average


%% Load probabilistic ROI (flat patch) to use for all subjs
ROI_lh = read_patch([fs_dir hemis{1} '.' latmed '_VisAudWM_combined_TFCE_flat.patch']);
ROI_rh = read_patch([fs_dir hemis{2} '.' latmed '_VisAudWM_combined_TFCE_flat.patch']);
ROIs = {ROI_lh, ROI_rh};

%% Compute group level axis of greatest change
group_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/grouplevel/';
new_x = {nan(2,1), nan(2,1)};
new_y = {nan(2,1), nan(2,1)};
group_data_diff = {};

for hh = 1:length(hemis)
    hemi = hemis{hh};

    % Precompute 100x100 mesh grid for plotting
    [xq,yq] = meshgrid(linspace(min(ROIs{hh}.x), max(ROIs{hh}.x), 100),...
                   linspace(min(ROIs{hh}.y), max(ROIs{hh}.y), 100));

    % Get group contrast data
    for cc = 1:length(contrasts)
        group_data = MRIread([group_dir hemis{hh} '.ces.localizer_groupavg_' contrasts{cc} '.glmres/osgm/z.mgh']);
        group_data_ROI{cc} = group_data.vol(ROIs{hh}.ind+1);
    end

    % Take difference of 2 contrasts
    group_data_diff{hh} = group_data_ROI{2} - group_data_ROI{1};
    
    switch axis_method
        case 'regression'
            % Linear regression to find 2D plane that best fits 3D data
            X = [ROIs{hh}.x; ROIs{hh}.y; ones(1,length(ROIs{hh}.x))]';
            coefs = X \ group_data_diff{hh}';
            coefs_norm = coefs(1:2) / norm(coefs(1:2)); % normalize coefficients
        
            % Plot 3D data surface
            zq = griddata(ROIs{hh}.x, ROIs{hh}.y, group_data_diff{hh}, xq, yq, 'linear'); % interpolate to 100x100 grid
            figure;
            surf(xq, yq, zq); % plot surface
            hold on;
            
            % Plot 2D plane of best fit
            zgrid = coefs(1) * xq + coefs(2) * yq + coefs(3); % Get z coord for each xy pair
            mesh(xq,yq,zgrid,'FaceColor', 'g', 'FaceAlpha',0.5, 'EdgeColor','none'); % plot plane
        
            % Plot 2D line along which there is the greatest change in z
            x0 = mean(ROIs{hh}.x);
            y0 = mean(ROIs{hh}.y);
            z0 = coefs(1)*x0 + coefs(2)*y0 + coefs(3);
            x_pts = linspace(min(ROIs{hh}.x)-x0, max(ROIs{hh}.x)-x0, 100);
            y_pts = linspace(min(ROIs{hh}.y)-y0, max(ROIs{hh}.y)-y0, 100);
            x_line = x0 + x_pts * coefs_norm(1);
            y_line = y0 + y_pts * coefs_norm(2);
            plot(x_line, y_line, 'r--', 'LineWidth',2);
            xlabel('x');
            ylabel('y');
            zlabel('PSC Difference');
        
            % Get new x and y axes directions where x axis is now 2d line of greatest change in z
            new_x{hh} = coefs_norm;
            new_y{hh} = [-coefs_norm(2); coefs_norm(1)]; % perpendicular to new_x (rotate 90 deg)
        case 'average'
            
    end
end

%% Rotate lh and rh flatmap coordinates to new x and y axes
for hh = 1:length(hemis)
    xy = [ROIs{hh}.x; ROIs{hh}.y];
    R = [new_x{hh}, new_y{hh}];
    xy_rotated = R * xy;
    ROIs{hh}.x = xy_rotated(1,:);
    ROIs{hh}.y = xy_rotated(2,:);
    
    % Plot psc data along new axes
    [xq,yq] = meshgrid(linspace(min(ROIs{hh}.x), max(ROIs{hh}.x), 100),...
                   linspace(min(ROIs{hh}.y), max(ROIs{hh}.y), 100));
    zq = griddata(ROIs{hh}.x, ROIs{hh}.y, group_data_diff{hh}, xq, yq, 'linear'); % interpolate to 100x100 grid
    figure;
    surf(xq, yq, zq); % plot surface    
    xlabel('x'); ylabel('y'); zlabel('PSC');
end


%% Loop over subjects and get contrast data within ROI
ROI_pscs = {nan(length(ROIs{1}.ind),2,N_subjs), nan(length(ROIs{2}.ind),2,N_subjs)};

for hh = 1:length(hemis)
    hemi = hemis{hh};
    [xq,yq] = meshgrid(linspace(min(ROIs{hh}.x), max(ROIs{hh}.x), 100),...
                   linspace(min(ROIs{hh}.y), max(ROIs{hh}.y), 100));

    % Calculate 10 possible boundary/gradient x-axis 'midpoints'
    x_step = range(ROIs{hh}.x)/11;
    midpts = min(ROIs{hh}.x)+x_step:x_step:max(ROIs{hh}.x)-x_step;

    % Make baseline boundary and gradient model for each midpoint
    boundary_model = {};
    gradient_model = {};
    for mm = 1:length(midpts)
        
    end

    for ss = 1:N_subjs
        subjCode = subjCodes{ss};
        for cc = 1:length(contrasts)
            % Get contrast PSC data
            contrast = contrasts{cc};
            psc = MRIread([subjdir subjCode '/localizer/localizer_contrasts_0sm_' hemi '/' contrast '/cespct.nii.gz']);
            ROI_psc = psc.vol(ROIs{hh}.ind+1);
            ROI_pscs{hh}(:,cc,ss) = ROI_psc;
        end
        
        % Get difference between contrasts in ROI
        ROI_psc_diffs = ROI_pscs{hh}(:,2,ss) - ROI_pscs{hh}(:,1,ss);
        
        % Loop through midpoints and fit each model
        for mm = 1:length(midpts)
            % Fit boundary model
        

            % Fit gradient model
        end

    end
end



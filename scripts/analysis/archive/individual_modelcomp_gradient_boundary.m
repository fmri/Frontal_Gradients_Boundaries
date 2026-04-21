%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to compare a gradient-based and
%%% boundary-based model of change in function (contrast t-stats) across an ROI in
%%% individual subjects
%%%
%%%
%%% Tom Possidente - Septemeber 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% NOTES %%%
% This 2D method may misclassify jagged edged boundaries as gradients due
% to ignoring one axis. To rectify this we could run these models in slices
% across the Y axis, and average out which model fit is the best across all
% slices. That way, if there is a jagged boundary, it will contriubte much
% less to variability in change of z across x axis. 
%
%%% NOTES %%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables
hemis = {'lh', 'rh'};
num_hemis = length(hemis);
fs_num = 163842; % number of vertices in fsaverage

% Get participant IDs
experiment_name = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode(~ismember(subjDf_cut.subjCode, {'AH', 'SL', 'RR', 'PP', 'MM'})); % rejected subjs (movement, no fixation, poor activation)
N_subjs = length(subjCodes);

% Set up directories
subjdir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';
label_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/probabilistic_allROIs/';
group_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/grouplevel/';

% Set models/methods used
contrasts = {'f-vP', 'vA-vP'}; % which functional data contrasts to use
keyword = 'vis'; % vis or aud
ROI_name = 'inf_lat_frontal'; % Which ROI to look at: midIFS, aINS, preSMA, inf_lat_frontal, sup_lat_frontal, cIPS
models = {'step', 'linear', 'hinge'};
axis_method = 'average'; % regression (uses regression to find axis of largest difference) or average (uses weighted average of positive and negative pts to make line)
axis_choice = 'step'; % step (use step function axis for all models) or individual (use best axis for each model individually)

plot_fits = true; % plot out individual subj fits

%% Load probabilistic ROI (flat patch and labels) to use for all subjs
% Label files used to map flat patch vertices back to RAS coordinates
label_lh = readtable([label_dir hemis{1} '.' ROI_name '_prob_thresh5.label'], 'FileType','text');
label_rh = readtable([label_dir hemis{2} '.' ROI_name '_prob_thresh5.label'], 'FileType','text');
labels = {label_lh, label_rh};

% Flat patch used in model fitting as x and y coords
ROI_lh = read_patch([label_dir hemis{1} '.' ROI_name '_prob_thresh5_flat.patch']);
ROI_rh = read_patch([label_dir hemis{2} '.' ROI_name '_prob_thresh5_flat.patch']);
ROIs = {ROI_lh, ROI_rh};

% Make sure all vertices in probabilistic patch and probabilistic label are the same (there can be small differences due to the way the patches are cut)
for hh = 1:2
    del_patch = ~ismember(ROIs{hh}.ind, labels{hh}{:,1});
    del_label = ~ismember(labels{hh}{:,1}, ROIs{hh}.ind);
    labels{hh}(del_label,:) = [];
    ROIs{hh}.ind = ROIs{hh}.ind(~del_patch);
    ROIs{hh}.x = ROIs{hh}.x(~del_patch);
    ROIs{hh}.y = ROIs{hh}.y(~del_patch);
    ROIs{hh}.z = ROIs{hh}.z(~del_patch);
    ROIs{hh}.vno = ROIs{hh}.vno(~del_patch);
    ROIs{hh}.npts = length(ROIs{hh}.ind);
end

%% Compute group level axis of greatest change
x_binedges = {nan(2,1), nan(2,1)};
new_y = {nan(2,1), nan(2,1)};
group_data_diff = {};

for hh = 1:num_hemis
    hemi = hemis{hh};

    % Precompute 100x100 mesh grid for plotting
    [xq,yq] = meshgrid(linspace(min(ROIs{hh}.x), max(ROIs{hh}.x), 100),...
        linspace(min(ROIs{hh}.y), max(ROIs{hh}.y), 100));

    % Get group contrast data
    for cc = 1:length(contrasts)
        group_data_curr = MRIread([group_dir hemis{hh} '.ces.localizer_groupavg_' contrasts{cc} '.glmres/osgm/z.mgh']);
        if strcmp(contrasts{cc}(1:2), 'f-') % this contrast got coded backwards, so flip it
            group_data_curr.vol = -group_data_curr.vol;
        end
        group_data_ROI{cc} = group_data_curr.vol(ROIs{hh}.ind+1); % indices are one off so add one
    end

    % Lower bound stats at 0 (more interpretable when taking the difference of 2 contrasts)
    group_data_ROI{1}(group_data_ROI{1}<0) = 0;
    group_data_ROI{2}(group_data_ROI{2}<0) = 0;

    % Take difference of 2 contrasts
    group_data_diff{hh} = group_data_ROI{2} - group_data_ROI{1};

    switch axis_method % 2 possible methods for finding axis of greatest change
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
            x_pts = linspace(min(ROIs{hh}.x)-x0, max(ROIs{hh}.x)-x0, 100); % will have to do this differently, 100);
            y_pts = linspace(min(ROIs{hh}.y)-y0, max(ROIs{hh}.y)-y0, 100);
            x_line = x0 + x_pts * coefs_norm(1);
            y_line = y0 + y_pts * coefs_norm(2);
            plot(x_line, y_line, 'r--', 'LineWidth',2);
            xlabel('x');
            ylabel('y');
            zlabel('T-stat Difference');

            % Get new x and y axes directions where x axis is now 2d line of greatest change in z
            x_binedges{hh} = coefs_norm;
            new_y{hh} = [-coefs_norm(2); coefs_norm(1)]; % perpendicular to new_x (rotate 90 deg)

        case 'average'
            % Get weighted average of x,y coords with positive z
            pos_mask = group_data_diff{hh} >= 0;
            weighted_pos(1) = sum(ROIs{hh}.x(pos_mask) .* group_data_diff{hh}(pos_mask)) / sum(group_data_diff{hh}(pos_mask));
            weighted_pos(2) = sum(ROIs{hh}.y(pos_mask) .* group_data_diff{hh}(pos_mask)) / sum(group_data_diff{hh}(pos_mask));

            % Get weighted average of x,y coords with negative z
            neg_mask = group_data_diff{hh} < 0;
            weighted_neg(1) = sum(ROIs{hh}.x(neg_mask) .* group_data_diff{hh}(neg_mask)) / sum(group_data_diff{hh}(neg_mask));
            weighted_neg(2) = sum(ROIs{hh}.y(neg_mask) .* group_data_diff{hh}(neg_mask)) / sum(group_data_diff{hh}(neg_mask));

            % Direction vector from pos to neg means will be new x axis
            v = weighted_pos - weighted_neg;
            x_binedges{hh} = v' / norm(v);  % normalize

            % Perpendicular vector (90 deg rotation)
            new_y{hh} = [-x_binedges{hh}(2); x_binedges{hh}(1)];

            % Plot 3D data surface
            figure;
            scatter3(ROIs{hh}.x, ROIs{hh}.y, group_data_diff{hh});
            hold on;
            plot([weighted_pos(1), weighted_neg(1)], [weighted_pos(2), weighted_neg(2)], 'LineWidth',5, 'Color', 'r');
            xlabel('x'); ylabel('y'); zlabel('T-stats')
            title([hemis{hh} ' flatmap with axis plotted'])
    end
end

%% Rotate lh and rh group flatmap coordinates to new x and y axes of greatest change
x_axis_rotation = nan(2,1);
for hh = 1:num_hemis
    xy = [ROIs{hh}.x; ROIs{hh}.y];
    R = [x_binedges{hh}, new_y{hh}]'; % rotation matrix
    xy_rotated = R * xy;
    x_axis_rotation(hh) = rad2deg(atan2(R(2,1), R(1,1))); % store degree of rotation
    ROIs{hh}.x = xy_rotated(1,:); % new coords
    ROIs{hh}.y = xy_rotated(2,:); % new coords

    figure;
    scatter(ROIs{hh}.x, ROIs{hh}.y, [], group_data_diff{hh});
    clim_set = max(abs(prctile(group_data_diff{hh}, [10,90])));
    xlabel('x'); ylabel('y'); cb = colorbar; colormap(redbluedark); clim([-clim_set, clim_set]); ylabel(cb, 'T-stat', 'rotation', 270);
    title([hemis{hh} ' flatmap rotated x greatest change'])
end


%% Loop over subjects and fit boundary and gradient models to contrast data
% Lots of stats/diagnostic storage variables
ROI_Ts = {nan(length(ROIs{1}.ind),2,N_subjs), nan(length(ROIs{2}.ind),2,N_subjs)};
winning_model = nan(N_subjs, num_hemis);
model_comp = nan(N_subjs, num_hemis);
winning_rsquare = nan(N_subjs, num_hemis);
linear_xdist = nan(N_subjs, num_hemis);
BICs = nan(N_subjs, num_hemis, 3);
hinge_direction = nan(N_subjs, num_hemis);

% Parameters for how to rotate axis when fitting
deg_step = 1; % degrees
deg_tolerance = 45;
angles = 0-deg_tolerance:deg_step:0+deg_tolerance;
winning_angle = nan(N_subjs, num_hemis, 3);

for hh = 1:num_hemis
    hemi = hemis{hh};

    for ss = 1:N_subjs
        subjCode = subjCodes{ss};
        for cc = 1:length(contrasts)
            % Get contrast t-stat data
            contrast = contrasts{cc};
            if strcmp(contrast(1:2), 'f-') % this contrast got coded backwards, so flip it
                contrast = [contrast(3:end) '-f'];
            end
            tstats = MRIread([subjdir subjCode '/localizer/localizer_contrasts_0sm_' hemi '/' contrast '/t.nii.gz']);
            ROI_tstat = tstats.vol(ROIs{hh}.ind+1); % indices are one off so add one
            ROI_Ts{hh}(:,cc,ss) = ROI_tstat;
        end

        % Clip negative t-stats at zero so that subtraction of 2 contrasts makes sense
        ROI_Ts{hh}(ROI_Ts{hh}(:,1,ss)<0,1,ss) = 0;
        contrast1_Ts = ROI_Ts{hh}(:,1,ss);
        ROI_Ts{hh}(ROI_Ts{hh}(:,2,ss)<0,2,ss) = 0;
        contrast2_Ts = ROI_Ts{hh}(:,2,ss);

        % Get difference between contrasts in ROI
        ROI_t_diffs = ROI_Ts{hh}(:,2,ss)  - ROI_Ts{hh}(:,1,ss);

        if plot_fits
            figure; subplot(2,1,1);
            scatter(ROIs{hh}.x, ROIs{hh}.y, [], ROI_t_diffs);
            clim_set = max(abs(prctile(ROI_t_diffs, [10,90])));
            xlabel('x'); ylabel('y'); cb = colorbar; colormap(redbluedark); clim([-clim_set, clim_set]); ylabel(cb, 'T-stat', 'rotation', 270);
            xlim([min(ROIs{hh}.x),max(ROIs{hh}.x)]); ylim([min(ROIs{hh}.y),max(ROIs{hh}.y)]);
        end

        % Replace outliers with 3 stds away from median value (clipping)
        ROI_t_diffs = filloutliers(ROI_t_diffs, 'clip', 'median'); % 3 stds from median is outlier

        if plot_fits
            subplot(2,1,2);
            scatter(ROIs{hh}.x, ROIs{hh}.y, [], ROI_t_diffs);
            clim_set = max(abs(prctile(ROI_t_diffs, [10,90])));
            xlabel('x'); ylabel('y'); cb = colorbar; colormap(redbluedark); clim([-clim_set, clim_set]); ylabel(cb, 'T-stat', 'rotation', 270);
            xlim([min(ROIs{hh}.x),max(ROIs{hh}.x)]); ylim([min(ROIs{hh}.y),max(ROIs{hh}.y)]);
        end
        
        % Rename for clarity/ease
        xs = ROIs{hh}.x;
        ys = ROIs{hh}.y;
        Ts = ROI_t_diffs;
        coord_inds = ROIs{hh}.ind; % freesurfer indices of these vertices (so we can map them back to RAS eventually)
        group_data_curr = group_data_diff{hh};

        % Start plot
        if plot_fits
            f1 = figure;
            subplot(2,2,1);
            scatter(xs, Ts);
            title('original');
        end

        %% Fit step-wise funciton to data (boundary model)
        [gof_step, info_step, winning_angle(ss,hh,1), xs_step, ys_step, fit_res_step] = fit_GB_model('step', xs, ys, Ts, group_data_curr, angles);

        %step_direction(ss,hh) = step_results_best(find(min(xs_step)==xs_step)) < step_results_best(find(max(xs_step)==xs_step));

        if plot_fits
            subplot(2,2,2);
            scatter3(xs_step, ys_step, Ts); xlabel('x'); ylabel('y'); zlabel('t-stats'); hold on; %scatter(xs_step, step_results_best);
            plot(fit_res_step)
            view([0 0]);
            title(['step model | rotate ' num2str(winning_angle(ss,hh,1)) ' | rsqr: ' num2str(round(gof_step.rsquare,3))])
        end

        %% Fit 2D plane to data (linear gradient model)
        
        if strcmp(axis_choice, 'step')
            angles_linear = winning_angle(ss,hh,1);
        else
            angles_linear = angles;
        end

        [gof_linear, info_linear, winning_angle(ss,hh,2), xs_linear, ys_linear, fit_res_linear] = ...
            fit_GB_model('linear', xs, ys, Ts, group_data_curr, angles_linear);

        if plot_fits
            subplot(2,2,3);
            scatter3(xs_linear, ys_linear, Ts); xlabel('x'); ylabel('y'); zlabel('t-stats'); hold on; %
            plot(fit_res_linear);
            view([0 0]);
            title(['linear model | rotate ' num2str(winning_angle(ss,hh,2)) ' | rsqr: ' num2str(round(gof_linear.rsquare,3))])
        end

        %% Fit 2D plane to subsection of data (hnige function gradient model)

        if strcmp(axis_choice, 'step')
            angles_hinge = winning_angle(ss,hh,1);
        else
            angles_hinge = angles;
        end

        [gof_hinge, info_hinge, winning_angle(ss,hh,3), xs_hinge, ys_hinge, fit_res_hinge] = ...
            fit_GB_model('hinge', xs, ys, Ts, group_data_curr, angles_hinge);

        hinge_direction(ss,hh) = fit_res_hinge.a < fit_res_hinge.b;

        if plot_fits
            subplot(2,2,4);
            scatter3(xs_hinge, ys_hinge, Ts); xlabel('x'); ylabel('y'); zlabel('t-stats'); hold on; %
            plot(fit_res_hinge);
            view([0 0]);
            title(['linear hinge | rotate ' num2str(winning_angle(ss,hh,3)) ' | rsqr: ' num2str(round(gof_hinge.rsquare,3))])
        end

        %% Calculate log likelihood of each model
        LL_step = loglikelihood(info_step.numobs, gof_step.rmse, gof_step.sse);
        LL_linear = loglikelihood(info_linear.numobs, gof_linear.rmse, gof_linear.sse);
        LL_hinge = loglikelihood(info_hinge.numobs, gof_hinge.rmse, gof_hinge.sse);

        % Compare AIC/BIC between models
        [aic, bic] = aicbic([LL_step; LL_linear; LL_hinge], [info_step.numparam; info_linear.numparam; info_hinge.numparam]);
        [~,ind] = min(bic); % get ind of min BIC model
        disp([models{ind}]) % print winning model name
        winning_model(ss,hh) = ind;
        model_comp(ss,hh) = bic(ind) - min(bic(~ismember(1:3,ind)));  % compare best model BIC to next best model BIC
        BICs(ss,hh,:) = bic; % record all BICs

        % Check if the winning model fit the data well
        gofs = {gof_step, gof_linear, gof_hinge};
        disp(['winning model r^2: ' num2str(round(gofs{ind}.rsquare,3))]);
        winning_rsquare(ss,hh) = gofs{ind}.rsquare;

        % If hinge model is the winner, check whether the linear piece is large enough to cross multiple voxels
        if ind==3
            mean_y1 = mean(ys_hinge( (xs_hinge<fit_res_hinge.x1+0.5) & (xs_hinge>fit_res_hinge.x1-0.5) ) ); % get mean y coord near x1
            mean_y2 = mean(ys_hinge( (xs_hinge<fit_res_hinge.x2+0.5) & (xs_hinge>fit_res_hinge.x2-0.5) ) ); % get mean y coord near x2
            linear_extent = [fit_res_hinge.x1, mean_y1, ; fit_res_hinge.x2, mean_y2]; % record coords at beginning and end of hinge

            % Find vertex nearest (x1, mean_y) and (x2, mean_y) and get RAS coords
            [dist1, vert1_ind] = min(cell2mat(arrayfun(@(x) pdist([linear_extent(1,:); [xs_hinge(x), ys_hinge(x)]]), 1:length(xs_hinge), 'UniformOutput', false) ) );
            if dist1>0.5
                disp('Poor vertex match');
                keyboard;
            end
            vert1 = coord_inds(vert1_ind); % which vertex index is nearest to point at beginning of hinge?
            RAS_vert1 = labels{hh}{labels{hh}.Var1==vert1, 2:4}; % find RAS coordinates of that vertex

            [dist2, vert2_ind] = min(cell2mat(arrayfun(@(x) pdist([linear_extent(2,:); [xs_hinge(x), ys_hinge(x)]]), 1:length(xs_hinge), 'UniformOutput', false) ) );
            if dist2>0.5
                disp('Poor vertex match');
                keyboard;
            end 
            vert2 = coord_inds(vert2_ind); % which vertex index is nearest to point at end of hinge?
            RAS_vert2 = labels{hh}{labels{hh}.Var1==vert2, 2:4}; % find RAS coordinates of that vertex

            % Find distance between vertices in RAS space
            linear_xdist(ss,hh) = pdist([RAS_vert1; RAS_vert2]);
        end

        if plot_fits
            sgtitle([subjCode ' ' hemi ' | Winner:' num2str(winning_model(ss,hh)) ' | BIC Diff:' num2str(round(model_comp(ss,hh),2)) ' | r^2:' num2str(round(winning_rsquare(ss,hh),3)) ' | hingedist: ' num2str(round(linear_xdist(ss,hh),2))])
        end

    end
end


%% Basic counts to compare models
sum(winning_model==3)
sum(winning_model==3 & hinge_direction==1)
sum(winning_model==3 & winning_rsquare>0.1)
sum(winning_model==3 & model_comp<-10)
sum(winning_model==3 & linear_xdist>4)

total_good = sum(winning_model==3 & winning_rsquare>0.1 & model_comp<-10 & linear_xdist>4 & hinge_direction==1)

mean(winning_rsquare(winning_model==3))
std(winning_rsquare(winning_model==3))
std(winning_rsquare(winning_model==3))/sqrt(sum(winning_model==3,'all'))

mean(linear_xdist(winning_model==3))
std(linear_xdist(winning_model==3))
std(linear_xdist(winning_model==3))/sqrt(sum(winning_model==3,'all'))

%% Group plots/tests
BIC_allhemis = reshape(BICs, [size(BICs,1)*size(BICs,2), size(BICs,3)]);
hingedir_allhemis = reshape(hinge_direction, [size(hinge_direction,1)*size(hinge_direction,2),1]);
winning_model_allhemis = reshape(winning_model, [size(winning_model,1)*size(winning_model,2),1]);

% Plot BIC differences
new_N = N_subjs*2-sum(hingedir_allhemis==0);
m = mean(BIC_allhemis(hingedir_allhemis==1,1)-BIC_allhemis(hingedir_allhemis==1,3));
figure; 
subplot(1,3,1);
swarmchart(ones(new_N,1), BIC_allhemis(hingedir_allhemis==1,1)-BIC_allhemis(hingedir_allhemis==1,3), [], 'c', 'filled'); hold on;
swarmchart(1, m, 100, 'r', 'filled');
yline(10, '--r');
xlim([-2,4])
ylim([-100,2500]);
ylabel('BIC Difference');
xticks([]);
title({['Step vs Hinge Model | \mu=' num2str(round(m,2))]})

% Plot hinge distances
xdists_hingewinners = linear_xdist(winning_model==3 & hinge_direction==1);
[h,p,CI,stats] = ttest(xdists_hingewinners(:), 4);
subplot(1,3,2);
swarmchart(ones(length(xdists_hingewinners)), xdists_hingewinners, [], 'c', 'filled'); hold on;
swarmchart(1, mean(xdists_hingewinners), 100, 'r', 'filled');
yline(4, '--r');
xlim([-2,4])
ylabel('Hinge Length Across Cortex (mm)');
xticks([]);
title({['Hinge Length | \mu=' num2str(round(mean(xdists_hingewinners),2))]})

% Plot hinge R squareds
rsqr_hingewinners = winning_rsquare(winning_model==3 & hinge_direction==1);
subplot(1,3,3);
swarmchart(ones(length(rsqr_hingewinners)), rsqr_hingewinners, [], 'c', 'filled'); hold on;
swarmchart(1, mean(rsqr_hingewinners), 100, 'r', 'filled');
yline(0.1, '--r');
xlim([-2,4])
ylabel('R-squared');
xticks([]);
title(['Hinge Model R-squared | \mu=', num2str(round(mean(rsqr_hingewinners),2))])

sgtitle([ROI_name ' ' keyword ' | N = ' num2str(new_N) ', ' num2str(length(rsqr_hingewinners))]);

%% Helper functions %%
function LL= loglikelihood(num_obs, rmse, sse)
    LL = -0.5 * num_obs * log( 2*pi*(rmse^2) ) - (1/(2*(rmse^2))) * sse;
end




%% Code Graveyard
% Limit to areas with significant signal in x and y axes (from sides, but not in middle)
% if ~strcmp(ROI_name, 'PVC')
%     binsize = 3; %mm
%     t_thresh = 1;
%     xbad_signal = limit_ROI_edges(ROIs{hh}.x, contrast1_Ts, contrast2_Ts, binsize, t_thresh);
%     ybad_signal = limit_ROI_edges(ROIs{hh}.y, contrast1_Ts, contrast2_Ts, binsize, t_thresh);
%     bad_signal = xbad_signal | ybad_signal;
% else
%     bad_signal = false(1, length(contrast1_Ts));
% end

% Nan out regions on tails of x and y axes that are not strong signal
% ROI_t_diffs(bad_signal) = nan;

% Get limited x,y coords and values for current ROI
% xs = ROIs{hh}.x(~bad_signal);
% ys = ROIs{hh}.y(~bad_signal);
% Ts = ROI_t_diffs(~bad_signal);
% coord_inds = ROIs{hh}.ind(~bad_signal); % freesurfer indices of these vertices (so we can map them back to RAS eventually)

% Cut group data to the same indices in order to derive parameter guesses for models from group data
% group_data_curr = group_data_diff{hh}(~bad_signal);

% function bad_inds = limit_ROI_edges(coords, stats1, stats2, binsize, stat_thresh)
% 
%     bin_edges = min(coords):binsize:max(coords);
%     nbins = length(bin_edges)-1;
%     binvec = true(nbins,1);
%     for bb = 1:nbins
%         inds = coords>bin_edges(bb) & coords<=bin_edges(bb+1);
%         data_inbin1 = mean(stats1(inds));
%         data_inbin2 = mean(stats2(inds));
%         if ~(data_inbin1 >= stat_thresh || data_inbin2 >= stat_thresh) % if neither bin has T>thresh mean, deselect that area
%             binvec(bb) = false;
%         end
%     end
% 
%     % Find start and stop inds for bins
%     start_ind = 1;
%     if binvec(1)~=1
%         start_ind = find(diff(binvec)==1, 1, 'first') + 1;
%     end
%     end_ind = nbins;
%     if binvec(nbins)~=1
%         end_ind = find(diff(binvec)==-1, 1, 'last');
%     end
% 
%     bad_inds = coords<bin_edges(start_ind) | coords>bin_edges(end_ind+1);
% end

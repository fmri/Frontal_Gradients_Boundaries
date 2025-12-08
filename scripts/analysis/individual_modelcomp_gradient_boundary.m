%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to compare a gradient-based and
%%% boundary-based model of change in function (contrast t-stats) across an ROI in
%%% individual subjects
%%%
%%%
%%% Tom Possidente - Septemeber 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
subjCodes = subjDf_cut.subjCode(~ismember(subjDf_cut.subjCode, {'AH', 'SL', 'RR', 'PP', 'MM'}));
N_subjs = length(subjCodes);

% Set up directories
subjdir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';
fs_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/fsaverage/surf/';
label_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/';
group_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/grouplevel/';

% Set models/methods used
contrasts = {'aPvP-f', 'vAaA-vPaP'}; % which functional data contrasts to use
%contrasts = {'f-vP', 'vA-vP'}; % which functional data contrasts to use
ROI_name = 'medial'; % Which ROI to look at: PVC, lateral_middle, lateral_inferior, lateral_superior, or medial
models = {'step', 'linear', 'hinge'};
axis_method = 'average'; % regression (uses regression to find axis of largest difference) or average (uses weighted average of positive and negative pts to make line)

plot_fits = false; % plot out individual subj fits

%% Load probabilistic ROI (flat patch and labels) to use for all subjs
% Label files used to map flat patch vertices back to RAS coordinates
label_lh = readtable([label_dir hemis{1} '.' ROI_name '_VisAudWM_combined_TFCE.label'], 'FileType','text');
label_rh = readtable([label_dir hemis{2} '.' ROI_name '_VisAudWM_combined_TFCE.label'], 'FileType','text');
labels = {label_lh, label_rh};

% Flat patch used in model fitting as x and y coords
ROI_lh = read_patch([fs_dir hemis{1} '.' ROI_name '_VisAudWM_combined_TFCE_flat.patch']);
ROI_rh = read_patch([fs_dir hemis{2} '.' ROI_name '_VisAudWM_combined_TFCE_flat.patch']);
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
        group_data_cut = MRIread([group_dir hemis{hh} '.ces.localizer_groupavg_' contrasts{cc} '.glmres/osgm/z.mgh']);
        if strcmp(contrasts{cc}(1:2), 'f-') % this contrast got coded backwards, so flip it
            group_data_cut.vol = -group_data_cut.vol;
        end
        group_data_ROI{cc} = group_data_cut.vol(ROIs{hh}.ind+1); % indices are one off so add one
    end

    % Lower bound stats at 0 (more interpretable when taking the difference of 2 contrasts)
    group_data_ROI{1}(group_data_ROI{1}<0) = 0;
    group_data_ROI{2}(group_data_ROI{2}<0) = 0;

    % Take difference of 2 contrasts
    group_data_diff{hh} = group_data_ROI{2} - group_data_ROI{1};

    % Plot data in RAS coordinates
    figure;
    scatter3(labels{hh}{:,2}, labels{hh}{:,3}, labels{hh}{:,4});
    xlabel('R'); ylabel('A'); zlabel('S')
    title([hemis{hh} ' RAS coordinates'])

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
            x_pts = linspace(min(ROIs{hh}.x)-x0, max(ROIs{hh}.x)-x0, 100);
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
    R = [x_binedges{hh}, new_y{hh}]';
    xy_rotated = R * xy;
    x_axis_rotation(hh) = rad2deg(atan2(R(2,1), R(1,1)));
    ROIs{hh}.x = xy_rotated(1,:);
    ROIs{hh}.y = xy_rotated(2,:);

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
deg_tolerance = 30;
angles = 0-deg_tolerance:deg_step:0+deg_tolerance;
num_angles = length(angles);
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

        % Limit to areas with significant signal in x and y axes (from sides, but not in middle)
        if ~strcmp(ROI_name, 'PVC')
            binsize = 3; %mm
            t_thresh = 1;
            xbad_signal = limit_ROI_edges(ROIs{hh}.x, contrast1_Ts, contrast2_Ts, binsize, t_thresh);
            ybad_signal = limit_ROI_edges(ROIs{hh}.y, contrast1_Ts, contrast2_Ts, binsize, t_thresh);
            bad_signal = xbad_signal | ybad_signal;
        else
            bad_signal = false(1, length(contrast1_Ts));
        end

        % Get difference between contrasts in ROI
        ROI_t_diffs = ROI_Ts{hh}(:,2,ss)  - ROI_Ts{hh}(:,1,ss);

        if plot_fits
            figure; subplot(2,1,1);
            scatter(ROIs{hh}.x, ROIs{hh}.y, [], ROI_t_diffs);
            clim_set = max(abs(prctile(ROI_t_diffs, [10,90])));
            xlabel('x'); ylabel('y'); cb = colorbar; colormap(redbluedark); clim([-clim_set, clim_set]); ylabel(cb, 'T-stat', 'rotation', 270);
            xlim([min(ROIs{hh}.x),max(ROIs{hh}.x)]); ylim([min(ROIs{hh}.y),max(ROIs{hh}.y)]);
        end

        % Nan out regions on tails of x and y axes that are not strong signal
        ROI_t_diffs(bad_signal) = nan;

        % Replace outliers with 3 stds away from median value (clipping)
        ROI_t_diffs = filloutliers(ROI_t_diffs, 'clip', 'median'); % 3 stds from median is outlier

        if plot_fits
            subplot(2,1,2);
            scatter(ROIs{hh}.x, ROIs{hh}.y, [], ROI_t_diffs);
            clim_set = max(abs(prctile(ROI_t_diffs, [10,90])));
            xlabel('x'); ylabel('y'); cb = colorbar; colormap(redbluedark); clim([-clim_set, clim_set]); ylabel(cb, 'T-stat', 'rotation', 270);
            xlim([min(ROIs{hh}.x),max(ROIs{hh}.x)]); ylim([min(ROIs{hh}.y),max(ROIs{hh}.y)]);
        end

        % Get limited x,y coords and values for current ROI
        xs = ROIs{hh}.x(~bad_signal);
        ys = ROIs{hh}.y(~bad_signal);
        Ts = ROI_t_diffs(~bad_signal);
        coord_inds = ROIs{hh}.ind(~bad_signal); % freesurfer indices of these vertices (so we can map them back to RAS eventually)

        % Cut group data to the same indices in order to derive parameter guesses for models from group data
        group_data_cut = group_data_diff{hh}(~bad_signal);

        % Start plot
        if plot_fits
            f1 = figure;
            subplot(2,2,1);
            scatter(xs, Ts);
            title('original');
        end

        %% Fit step-wise funciton to data (boundary model)

        [gof_step, info_step, winning_angle(ss,hh,1), xs_step, ~, step_results_best] = fit_model('step', xs, ys, Ts, group_data_cut, angles);

        if plot_fits
            subplot(2,2,2);
            scatter(xs_step, Ts); hold on; scatter(xs_step, step_results_best);
            title(['step model | rotate ' num2str(winning_angle(ss,hh,1)) ' | rsqr: ' num2str(round(gof_step.rsquare,3))])
        end

        %% Fit 2D plane to data (gradient model)

        [gof_linear, info_linear, winning_angle(ss,hh,2), xs_linear, ys_linear, fit_res_linear] = ...
            fit_model('linear', xs, ys, Ts, group_data_cut, angles);

        if plot_fits
            subplot(2,2,3);
            scatter3(xs_linear, ys_linear, Ts); xlabel('x'); ylabel('y'); zlabel('t-stats'); hold on; %
            plot(fit_res_linear);
            view([0 0]);
            title(['linear model | rotate ' num2str(winning_angle(ss,hh,2)) ' | rsqr: ' num2str(round(gof_linear.rsquare,3))])
        end

        %% Fit 2D plane to subsection of data (gradient adjusted model)

        [gof_hinge, info_hinge, winning_angle(ss,hh,3), xs_hinge, ys_hinge, fit_res_hinge] = ...
            fit_model('hinge', xs, ys, Ts, group_data_cut, angles);


        hinge_direction(ss,hh) = fit_res_hinge.a < fit_res_hinge.b;

        if plot_fits
            subplot(2,2,4);
            scatter3(xs_hinge, ys_hinge, Ts); xlabel('x'); ylabel('y'); zlabel('t-stats'); hold on; %
            plot(fit_res_hinge);
            view([0 0]);
            title(['linear hinge | rotate ' num2str(winning_angle(ss,hh,3)) ' | rsqr: ' num2str(round(gof_hinge.rsquare,3))])
        end

        %% Calculate log likelihood of each model
        LL_step = -0.5 * info_step.numobs * log( 2*pi*(gof_step.rmse^2) ) - (1/(2*(gof_step.rmse^2))) * gof_step.sse;
        LL_linear = -0.5 * info_linear.numobs * log( 2*pi*(gof_linear.rmse^2) ) - (1/(2*(gof_linear.rmse^2))) * gof_linear.sse;
        LL_hinge = -0.5 * info_hinge.numobs * log( 2*pi*(gof_hinge.rmse^2) ) - (1/(2*(gof_hinge.rmse^2))) * gof_hinge.sse;

        % Compare AIC/BIC between models
        [aic, bic] = aicbic([LL_step; LL_linear; LL_hinge], [info_step.numparam; info_linear.numparam; info_hinge.numparam]);
        [~,ind] = min(bic); % get ind of min BIC model
        disp([models{ind}])
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
            mean_y2 = mean(ys_hinge( (xs_hinge<fit_res_hinge.x2+0.5) & (xs_hinge>fit_res_hinge.x2-0.5) ) ); % get mean y coord near x1
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


%% Plots to compare models

% Winning model bar graph
figure;
data_lh = winning_model(:,1);
data_rh = winning_model(:,2);
bar([1,2,3], [sum(data_lh==1), sum(data_lh==2), sum(data_lh==3); sum(data_rh==1), sum(data_rh==2), sum(data_rh==3)]);
xticklabels({'Step', 'Linear', 'Hinge'});
ylabel('Winning Model');

% Hinge model BIC differences
figure;
swarmchart(ones(sum(data_lh==3),1), model_comp(data_lh==3)); hold on;
swarmchart(ones(sum(data_rh==3),1)*2, model_comp(data_rh==3));
xticks([1,2]);
xticklabels({'lh', 'rh'});
ylabel('BIC Difference');
yline(-10,'--r');

% Hinge model r-squared
figure;
swarmchart(ones(sum(data_lh==3),1), winning_rsquare(data_lh==3)); hold on;
swarmchart(ones(sum(data_rh==3),1)*2, winning_rsquare(data_rh==3));
xticks([1,2]);
xticklabels({'lh', 'rh'});
ylabel('Model Fit (r-square)');
yline(0.05,'--r');

% RAS coord distance for winning hinge models
figure;
swarmchart(ones(sum(data_lh==3),1), linear_xdist(data_lh==3)); hold on;
swarmchart(ones(sum(data_rh==3),1)*2, linear_xdist(data_rh==3));
xticks([1,2]);
xticklabels({'lh', 'rh'});
ylabel('hinge distance (mm)');
yline(4,'--r');


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

%% Group tests
BIC_sums = squeeze(sum(BICs, [1,2]));
xdists_hingewinners = linear_xdist(winning_model==3);
[h,p,CI,stats] = ttest(xdists_hingewinners, 4);


%% Helper functions %%
function bad_inds = limit_ROI_edges(coords, stats1, stats2, binsize, stat_thresh)

    bin_edges = min(coords):binsize:max(coords);
    nbins = length(bin_edges)-1;
    binvec = true(nbins,1);
    for bb = 1:nbins
        inds = coords>bin_edges(bb) & coords<=bin_edges(bb+1);
        data_inbin1 = mean(stats1(inds));
        data_inbin2 = mean(stats2(inds));
        if ~(data_inbin1 >= stat_thresh || data_inbin2 >= stat_thresh) % if neither bin has T>thresh mean, deselect that area
            binvec(bb) = false;
        end
    end
    
    % Find start and stop inds for bins
    start_ind = 1;
    if binvec(1)~=1
        start_ind = find(diff(binvec)==1, 1, 'first') + 1;
    end
    end_ind = nbins;
    if binvec(nbins)~=1
        end_ind = find(diff(binvec)==-1, 1, 'last');
    end
    
    bad_inds = coords<bin_edges(start_ind) | coords>bin_edges(end_ind+1);
end

function [gof, info, winning_angle, xs_out, ys_out, model_out] = fit_model(type, xs, ys, Ts, group_data_cut, angles)

    best_rsquare = 0;
    
    if strcmp(type, 'step')
        step_model = @(params, x) params(2)*(x < params(1)) + params(3)*(x >= params(1)); % param(1) is step location, 2 is pre-step, 3 is post-step
    
        for aa = 1:length(angles) % rotate x-axis iteratively and fit each time to find best fit
            theta = deg2rad(angles(aa));
            xs_new = xs * cos(theta) - ys * sin(theta); % only have to rotate x axis for step model fit (does not use y coords)
            parameter_guesses = [mean(xs_new), mean(group_data_cut(group_data_cut<0)), mean(group_data_cut(group_data_cut<0))]; % estimates for step location, pre-step z value, and post-step z value
            step_fit = lsqcurvefit(step_model, parameter_guesses, xs_new, Ts'); % fit model to data
            step_results = step_model(step_fit, xs_new); % run fitted model
            residuals_step = Ts' - step_results; % calculate residuals
            sse = sum(residuals_step.^2);
            rsqr = 1 - (sse / sum( (Ts'-mean(Ts)).^2 ) );
            if best_rsquare<rsqr % best model has lowest r-square
                best_rsquare = rsqr;
                gof.rsquare = rsqr;
                winning_angle = angles(aa); % record winning angle
                gof.rmse = sqrt(mean(residuals_step.^2));
                gof.sse = sse;
                xs_out = xs_new; % record winning angle x coords
                model_out = step_results; % record winning angle results for plotting
            end
        end
        info.numobs = length(xs);
        info.numparam = 3;
        ys_out = []; % not needed for this model
    
        return
    
    elseif strcmp(type, 'linear')
        model = fittype('x*a + b + y*0',... % linear function
            'dependent', 'z',...
            'independent', {'x','y'}, ... % y is not actually used in the function but include it so that we can use all 3D data to fit (not just the x and z coordinates)
            'coefficients', {'a','b'}); % a is slope, b is intercept
    
        for aa = 1:length(angles) % rotate x-axis iteratively and fit each time to find best fit
            theta = deg2rad(angles(aa));
            xs_new = xs * cos(theta) - ys * sin(theta); % get new rotated x and y coords
            ys_new = xs * sin(theta) + ys * cos(theta);
    
            % Calculate initial guesses for gradient model parameters using group-level data
            X = [xs_new; ones(1,length(xs_new))]';
            coefs = X \ group_data_cut'; % linear regression on only x coordinates
            slope_guess = coefs(1); % 1st param is slope
            intercept_guess = coefs(2); % 2nd is intercept
            startpoint_guesses = [slope_guess, intercept_guess];
    
            [pre_fit_res, pre_gof, pre_info] = fit([xs_new', ys_new'], Ts, model, ...
                'StartPoint', startpoint_guesses);
    
            if best_rsquare<pre_gof.rsquare
                best_rsquare = pre_gof.rsquare;
                winning_angle = angles(aa);
                xs_out = xs_new;
                ys_out = ys_new;
                model_out = pre_fit_res;
                gof = pre_gof;
                info = pre_info;
            end
        end
    
    elseif strcmp(type, 'hinge')
        model = fittype('a*(x < x1) + b*(x > x2) + ( a + ( (b-a)/(x2-x1)  * (x-x1) ) ) * ( (x > x1) & (x < x2) ) + y*0', ... % piecewise function, constant from -Inf to x1, linear from x1 to x2, constant from x2 to Inf. Must include y in the function even if it has no mathematical effect
            'dependent', 'z',...
            'independent', {'x','y'}, ... % y is not actually used in the function but include it so that we can use all 3D data to fit (not just the x and z coordinates)
            'coefficients', {'a','b','x1', 'x2'}); % a is z from -Inf to x1, b is z from x2 to Inf, x1 is where to start linear piece, x2 is where to end linear piece
    
        for aa = 1:length(angles) % rotate x-axis iteratively and fit each time to find best fit
            theta = deg2rad(angles(aa));
            xs_new = xs * cos(theta) - ys * sin(theta); % get new rotated x and y coords
            ys_new = xs * sin(theta) + ys * cos(theta);
    
            % Calculate initial guesses for gradient model parameters using group-level data
            startpoint_guesses = [mean(group_data_cut(group_data_cut<0)), mean(group_data_cut(group_data_cut>0)), prctile(xs_new, 25), prctile(xs_new, 75)]; % estimates for
            lower_bounds = [-Inf, -Inf, min(xs_new), min(xs_new)];
            upper_bounds = [Inf, Inf, max(xs_new), max(xs_new)];
    
            [pre_fit_res, pre_gof, pre_info] = fit([xs_new', ys_new'], Ts, model, ...
                'StartPoint', startpoint_guesses,...
                'Lower', lower_bounds,...
                'Upper', upper_bounds); % fit model
    
            if pre_fit_res.x2< pre_fit_res.x1 % Invalid fit, retry with bounds that force x1<x2
                [pre_fit_res, pre_gof, pre_info] = fit([xs_new', ys_new'], Ts, model, ...
                    'StartPoint', startpoint_guesses,...
                    'Lower', [lower_bounds(1), lower_bounds(2), lower_bounds(3), pre_fit_res.x1],...
                    'Upper', [upper_bounds(1), upper_bounds(2), pre_fit_res.x2, upper_bounds(3)]);
                assert(pre_fit_res.x2>pre_fit_res.x1, 'invalid hinge fit');
            end
    
            if best_rsquare<pre_gof.rsquare
                best_rsquare = pre_gof.rsquare;
                winning_angle = angles(aa);
                xs_out = xs_new;
                ys_out = ys_new;
                model_out = pre_fit_res;
                gof = pre_gof;
                info = pre_info;
            end
        end
    
    else
        error('model type not recognized')
    end
end


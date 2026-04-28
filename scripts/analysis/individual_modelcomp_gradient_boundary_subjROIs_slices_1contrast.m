%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to compare a gradient-based and
%%% boundary-based model of change in function (contrast t-stats) across an ROI in
%%% individual subjects
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
subjCodes = {'MK','AB''AD','LA','AE','TP','NM','AF','AG','AI','GG','UV','PQ','KQ','LN','RT','PT','PL','NS'};
N_subjs = length(subjCodes);

% Set up directories
subjdir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';
label_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/probabilistic_allROIs/';
group_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/grouplevel/';
subj_ROIs = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/subj_specific_01/';

% Set models/methods used
contrast = 'aA-aP'; % which functional data contrasts to use 
modailty = 'auditory'; % visual, auditory, supramodal, visaud
ROI_name = 'inf_lat_frontal'; % Which ROI to look at: midIFS, aINS, preSMA, inf_lat_frontal, sup_lat_frontal, cIPS, visaud_cIFSG_midIFS_interface
models = {'step', 'linear', 'hinge'};
axis_method = 'regression'; % regression (uses regression to find axis of largest difference) or average (uses weighted average of positive and negative pts to make line)
axis_choice = 'step'; % step (use step function axis for all models) or individual (use best axis for each model individually)

plot_fits = false; % plot out individual subj fits (debugging only unless you want a ~100 plots)

dist_thresh = 4.4; % length in mm of 2 voxels

%% Load probabilistic ROI (flat patch) to use for determining group level axis of greatest change
ROI_lh = read_patch([label_dir hemis{1} '.' ROI_name '_prob_thresh5_flat.patch']);
ROI_rh = read_patch([label_dir hemis{2} '.' ROI_name '_prob_thresh5_flat.patch']);
ROIs = {ROI_lh, ROI_rh};

%% Compute group level axis of greatest change
group_data = cell(2,1);

for hh = 1:num_hemis
    hemi = hemis{hh};

    % Get group contrast data
    if strcmp(contrast(3:4), '-f') % the fixation contrast got coded backwards here so the data files are all 'f-something' instead of 'something-f'
        new_contrast_str = ['f-' contrast(1:2)]; % flip the string around
        group_data_curr = MRIread([group_dir hemis{hh} '.ces.localizer_groupavg_' new_contrast_str '.glmres/osgm/z.mgh']);
        group_data_curr.vol = -group_data_curr.vol; % flip data to get back to 'something-fixation'
    else
        group_data_curr = MRIread([group_dir hemis{hh} '.ces.localizer_groupavg_' contrast '.glmres/osgm/z.mgh']);
    end
    group_data_ROI = group_data_curr.vol(ROIs{hh}.ind+1); % index numbers in labels/patches start at 0, but indexing starts at 1 in matlab, so add one when indexing like this

    % Lower bound stats at 0 (more interpretable when taking the difference of 2 contrasts)
    group_data_ROI(group_data_ROI<0) = 0;
    group_data{hh} = group_data_ROI; 

    % Find axis of greatest change and rotate coordinates so X axis is greatest change axis
    [ROIs{hh}.x, ROIs{hh}.y] = group_find_axis(axis_method, group_data{hh}, ROIs{hh}.x', ROIs{hh}.y', hemi);
end

%% Loop over subjects and fit boundary and gradient models to contrast data
% Lots of stats/diagnostic storage variables
ROI_Ts = cell(N_subjs, 2); % T stats for each ROI/subj/hemi
winning_model = nan(N_subjs, num_hemis); % which model had best BIC (1 is step, 2 is linear, 3 is hinge)
model_comp = nan(N_subjs, num_hemis); % BIC difference values between winner and next best
winning_rsquare = nan(N_subjs, num_hemis); % r-square of winning model
linear_xdist = nan(N_subjs, num_hemis); % distance of hinge in hinge model (mm)
BICs = nan(N_subjs, num_hemis, 3); % Raw BICs for each model
change_direction = nan(N_subjs, num_hemis, 3); % 1 for expected direction (increase along X axis), else 0

% Parameters for how to rotate axis when fitting
deg_step = 1; % degrees
deg_tolerance = 60; % searching +/- 45 degrees off of group axis of greatest change
angles = 0-deg_tolerance:deg_step:0+deg_tolerance;
winning_angles = nan(N_subjs, num_hemis); % which angle was optimal

tic;
for hh = 1:num_hemis
    hemi = hemis{hh};

    for ss = 1:N_subjs
        subjCode = subjCodes{ss};

        % Get subj ROI
        subjROI_path = [subj_ROIs hemi '.' ROI_name '_' modailty '_' subjCode '.label'];
        if ~isfile(subjROI_path) % if this ROI doesn't exist, the subj didn't have enough good signal/vertices in this ROI, skip it
            continue;
        end

        subjROI = readtable(subjROI_path, 'FileType','text'); % load label file for individual subj ROI
        patch_inds_mask = ismember(ROIs{hh}.ind, subjROI{:,1}); % Get flat patch indices/x/y for each vertex in subj ROI
        patch_inds = ROIs{hh}.ind(patch_inds_mask);
        subjROI_x = ROIs{hh}.x(patch_inds_mask);
        subjROI_y = ROIs{hh}.y(patch_inds_mask);

        tstats = MRIread([subjdir subjCode '/localizer/localizer_contrasts_0sm_' hemi '/' contrast '/t.nii.gz']);
        ROI_tstat = tstats.vol(patch_inds+1); % index numbers in labels/patches start at 0, but indexing starts at 1 in matlab, so add one when indexing like this
        ROI_tstat(ROI_tstat<0) = 0; % Clip negative t-stats at zero so that subtraction of 2 contrasts makes sense
        ROI_Ts{ss,hh} = ROI_tstat;

        Ts = ROI_Ts{ss,hh};
        
        % Clip big outliers
        thresh = mean(Ts(Ts>0)) * 5;
        Ts(Ts > thresh) = thresh; 

        % Rename for clarity/ease
        xs = subjROI_x;
        ys = subjROI_y;
        group_data_curr = group_data{hh}(patch_inds_mask); % limit group data to inside subj ROI

        % Start plot
        if plot_fits % plot x-coords and T vals
            f1 = figure;
            subplot(2,2,1);
            scatter(xs, ys, [], Ts, 'filled');
            clim_set = max(abs(prctile(Ts, [10,90])));
            xlabel('x'); ylabel('y'); cb = colorbar; colormap('autumn'); clim([0, clim_set]); ylabel(cb, 'T-stat', 'rotation', 270);
            xlim([min(xs),max(xs)]); ylim([min(ys),max(ys)]);
            title('original');
        end

        %% Fit step-wise function to data (boundary model)

        winning_angle = individual_find_axis_noslice('step', xs, ys, Ts, group_data_curr, angles);
        winning_angles(ss,hh) = winning_angle; 
        [metrics_step, ~, ~, xs_new, ys_new, ~, bins_step] = fit_GB_model_slices('step', xs, ys, Ts, group_data_curr, winning_angle);

        if plot_fits % plot step model fit
            figure(f1);
            subplot(2,2,2);
            scatter(xs_new, ys_new, [], Ts, 'filled');
            clim_set = max(abs(prctile(Ts, [10,90])));
            xlabel('x'); ylabel('y'); cb = colorbar; colormap(redbluedark); clim([-clim_set, clim_set]); ylabel(cb, 'T-stat', 'rotation', 270);
            xlim([min(xs_new),max(xs_new)]); ylim([min(ys_new),max(ys_new)]);
            hold on;
            for bb = 1:length(bins_step)-1
                midbin = (bins_step(bb+1) + bins_step(bb)) / 2;
                yline(midbin, ':');
                scatter(metrics_step.p1(bb), midbin, 100, 'g', '^', 'filled')
            end
            title(['step: rotate= ' num2str(winning_angle) ', rsqr=' num2str(round(mean(metrics_step.rsquare, 'omitnan', "Weights", metrics_step.numobs),2))])
        end

        %% Fit 2D plane to data (linear gradient model)
        [metrics_linear, ~, ~, xs_new, ys_new, ~] = fit_GB_model_slices('linear', xs, ys, Ts, group_data_curr, winning_angle);

        if plot_fits % plot
            figure(f1);
            subplot(2,2,3);
            scatter(xs_new, ys_new, [], Ts, 'filled');
            clim_set = max(abs(prctile(Ts, [10,90])));
            xlabel('x'); ylabel('y'); cb = colorbar; colormap(redbluedark); clim([-clim_set, clim_set]); ylabel(cb, 'T-stat', 'rotation', 270);
            xlim([min(xs_new),max(xs_new)]); ylim([min(ys_new),max(ys_new)]);
            title(['linear: rsqr=' num2str(round(mean(metrics_linear.rsquare, 'omitnan', "Weights", metrics_step.numobs),2))])        
        end

        %% Fit 2D plane to subsection of data (hinge function gradient model)
        [metrics_hinge, ~, ~, xs_new, ys_new, ~] = fit_GB_model_slices('hinge', xs, ys, Ts, group_data_curr, winning_angle);

        if plot_fits % plot
            figure(f1);
            subplot(2,2,4);
            scatter(xs_new, ys_new, [], Ts, 'filled');
            clim_set = max(abs(prctile(Ts, [10,90])));
            xlabel('x'); ylabel('y'); cb = colorbar; colormap(redbluedark); clim([-clim_set, clim_set]); ylabel(cb, 'T-stat', 'rotation', 270);
            xlim([min(xs_new),max(xs_new)]); ylim([min(ys_new),max(ys_new)]);
            hold on;
            for bb = 1:length(bins_step)-1
                midbin = (bins_step(bb+1) + bins_step(bb)) / 2;
                yline(midbin, ':');
                scatter(metrics_hinge.p1(bb), midbin, 100, 'g', '^', 'filled')
                scatter(metrics_hinge.p1(bb)+metrics_hinge.p2(bb), midbin, 100, 'g', '>', 'filled')
            end
            title(['hinge: rsqr=' num2str(round(mean(metrics_hinge.rsquare, 'omitnan', "Weights", metrics_step.numobs),2))])        
        end
        %% Take weighted means of metrics
        step_means = mean(metrics_step, 'omitnan', 'Weights', metrics_step.numobs);
        linear_means = mean(metrics_linear,'omitnan', 'Weights', metrics_step.numobs);
        hinge_means = mean(metrics_hinge, 'omitnan', 'Weights', metrics_step.numobs);

        %% Calculate log likelihood of each model
        LL_step = loglikelihood(step_means.numobs, step_means.rmse, step_means.sse);
        LL_linear = loglikelihood(linear_means.numobs, linear_means.rmse, linear_means.sse);
        LL_hinge = loglikelihood(hinge_means.numobs, hinge_means.rmse, hinge_means.sse);

        % Compare BIC between models
        [~, bic] = aicbic([LL_step; LL_linear; LL_hinge], [3; 2; 4]);
        [~,ind] = min(bic); % get ind of min BIC model
        disp([models{ind}]) % print winning model name
        winning_model(ss,hh) = ind;
        model_comp(ss,hh) = bic(ind) - min(bic(~ismember(1:3,ind)));  % compare best model BIC to next best model BIC
        BICs(ss,hh,:) = bic; % record all BICs

        % Check if the winning model actually fit the data well
        rsquares = [step_means.rsquare, linear_means.rsquare, hinge_means.rsquare];
        disp(['winning model r^2: ' num2str(round(rsquares(ind),3))]);
        winning_rsquare(ss,hh) = rsquares(ind);
        if isinf(rsquares(ind)) || rsquares(ind) < 0
            keyboard;
        end

        % If hinge model is the winner, check whether the linear piece is large enough to cross multiple voxels
        if ind==3
            linear_xdist(ss,hh) = mean((metrics_hinge.p2 + metrics_hinge.p1) - metrics_hinge.p1, 'omitnan', 'Weights', metrics_hinge.numobs); % just take difference between end of hinge and beginning of hinge on x axis
        end

        if plot_fits % make overall title for 2x2 fit plot
            sgtitle([subjCode ' ' hemi ' | Winner:' num2str(winning_model(ss,hh)) ' | BIC Diff:' num2str(round(model_comp(ss,hh),2)) ' | r^2:' num2str(round(winning_rsquare(ss,hh),3)) ' | hingedist: ' num2str(round(linear_xdist(ss,hh),2))])
        end

        %% Check direction of change for each model
        change_direction(ss,hh,1) = mode(metrics_step.p2 < metrics_step.p3);
        change_direction(ss,hh,2) = mode(metrics_linear.p1 > 0);
        change_direction(ss,hh,3) = mode(metrics_hinge.p3 < metrics_hinge.p4);
    end
end
toc

% Check the change direction of each winning model
good_winning_direction = zeros(N_subjs,num_hemis);
for mm = 1:3
    good_winning_direction = good_winning_direction + ((winning_model==mm) .* change_direction(:,:,mm));
end

%% Basic counts to compare models
% Something is off here, these should all add up to total
total_gradient_win = sum(winning_model==3 & winning_rsquare>0.1 & model_comp<-10 & linear_xdist>dist_thresh & change_direction(:,:,1)==1) % put them all together
total_boundary_win = sum(( (winning_model==1 & model_comp<-10) | (winning_model==3 & linear_xdist<=dist_thresh) ) & winning_rsquare>0.1 & change_direction(:,:,3)==1) % put them all together
total_linear_win = sum(winning_model==2 & winning_rsquare>0.1 & model_comp<-10 & change_direction(:,:,2)==1) % put them all together
weak_evidence = sum( (~(winning_model==3 & linear_xdist<=dist_thresh) & model_comp>=-10) & winning_rsquare>0.1 & good_winning_direction==1)
bad_direction = sum(good_winning_direction==0 & winning_rsquare>0.1)
r2_rejects = sum(winning_rsquare<=0.1)

total = 36 - sum(isnan(winning_model), 'all')
total_check = sum([total_gradient_win,total_boundary_win,total_linear_win,weak_evidence,bad_direction,r2_rejects])

%% Group plots/tests
BIC_allhemis = reshape(BICs, [size(BICs,1)*size(BICs,2), size(BICs,3)]); % collapse hemispheres
hingedir_allhemis = reshape(change_direction(:,:,3), [size(change_direction,1)*size(change_direction,2),1]); % collapse hemispheres
winning_model_allhemis = reshape(winning_model, [size(winning_model,1)*size(winning_model,2),1]); % collapse hemispheres

% plot angles
figure;
histogram(winning_angles(:), 'NumBins',20);

% Plot BIC differences
new_N = (N_subjs*2) - sum(hingedir_allhemis==0 | isnan(hingedir_allhemis)); % calculate number of subjs/hemis with hinge in correct direction
BIC_diffs = BIC_allhemis(hingedir_allhemis==1,1)-BIC_allhemis(hingedir_allhemis==1,3); % when hinge was in correct direction, what was the BIC difference bt hinge and step?
m = mean(BIC_diffs);
SE = std(BIC_diffs)/sqrt(new_N);
figure;
subplot(1,3,1);
swarmchart(ones(new_N,1), BIC_diffs, [], 'k', 'filled'); hold on;
errorbar(1, m, SE, 'Color', 'r','Marker', '.', 'MarkerSize', 30, 'CapSize',15, 'LineWidth',3);
yline(10, '--r'); ylim([min(BIC_diffs)-20, max(BIC_diffs)+20]);
ylabel('BIC Difference');
xlim([-2,4]); xticks([]);
title({['Step vs Hinge Model | \mu=' num2str(round(m,2))]})

% Plot hinge distances
xdists_hingewinners = linear_xdist(winning_model==3 & change_direction(:,:,3)==1);
m = mean(xdists_hingewinners);
SE = std(xdists_hingewinners)/sqrt(length(xdists_hingewinners));
subplot(1,3,2);
swarmchart(ones(length(xdists_hingewinners)), xdists_hingewinners, [], 'k', 'filled'); hold on;
errorbar(1, m, SE, 'Color', 'r','Marker', '.', 'MarkerSize', 30, 'CapSize',15, 'LineWidth',3);
yline(dist_thresh, '--r');
xlim([-2,4]); xticks([]);
ylabel('Hinge Length Across Cortex (mm)');
title({['Hinge Length | \mu=' num2str(round(mean(xdists_hingewinners),2))]})

% Plot hinge R squareds
rsqr_hingewinners = winning_rsquare(winning_model==3 & change_direction(:,:,3)==1);
m = mean(rsqr_hingewinners);
SE = std(rsqr_hingewinners/sqrt(length(rsqr_hingewinners)));
subplot(1,3,3);
swarmchart(ones(length(rsqr_hingewinners)), rsqr_hingewinners, [], 'k', 'filled'); hold on;
errorbar(1, m, SE, 'Color', 'r','Marker', '.', 'MarkerSize', 30, 'CapSize',15, 'LineWidth',3);
yline(0.1, '--r');
xlim([-2,4]); xticks([]);
ylabel('R-squared');
title(['Hinge Model R-squared | \mu=', num2str(round(mean(rsqr_hingewinners),2))])

sgtitle([replace(ROI_name,'_','-') ' ' modailty ' | N = ' num2str(new_N) ', ' num2str(length(rsqr_hingewinners))]);

%% Helper functions %%
function LL= loglikelihood(num_obs, rmse, sse)
LL = -0.5 * num_obs * log( 2*pi*(rmse^2) ) - (1/(2*(rmse^2))) * sse;
end

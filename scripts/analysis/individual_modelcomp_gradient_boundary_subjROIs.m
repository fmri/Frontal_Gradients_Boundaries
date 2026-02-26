%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to compare a gradient-based and
%%% boundary-based model of change in function (contrast t-stats) across an ROI in
%%% individual subjects
%%%
%%% Tom Possidente - Septemeber 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% NOTES %%%
% This 2D method may misclassify jagged edged boundaries as gradients due
% to ignoring one axis. To rectify this we could run these models in slices
% across the Y axis, and average out which model fit is the best across all
% slices. That way, if there is a jagged boundary, it will contriubte much
% less to variability in change of z across x axis. 
%%% NOTES %%%

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
subj_ROIs = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/subj_specific_05/';

% Set models/methods used
contrasts = {'vP-f', 'vA-vP'}; % which functional data contrasts to use
N_contrasts = length(contrasts);
modailty = 'visual'; % vis or aud
ROI_name = 'preSMA'; % Which ROI to look at: midIFS, aINS, preSMA, inf_lat_frontal, sup_lat_frontal, cIPS
models = {'step', 'linear', 'hinge'};
axis_method = 'average'; % regression (uses regression to find axis of largest difference) or average (uses weighted average of positive and negative pts to make line)
axis_choice = 'step'; % step (use step function axis for all models) or individual (use best axis for each model individually)

plot_fits = true; % plot out individual subj fits (debugging only unless you want a ~100 plots)

%% Load probabilistic ROI (flat patch) to use for determining group level axis of greatest change
ROI_lh = read_patch([label_dir hemis{1} '.' ROI_name '_prob_thresh5_flat.patch']);
ROI_rh = read_patch([label_dir hemis{2} '.' ROI_name '_prob_thresh5_flat.patch']);
ROIs = {ROI_lh, ROI_rh};

%% Compute group level axis of greatest change
group_data_diff = cell(2,1);

for hh = 1:num_hemis
    hemi = hemis{hh};

    % Get group contrast data
    for cc = 1:length(contrasts)
        if strcmp(contrasts{cc}(3:4), '-f') % the fixation contrast got coded backwards here so the data files are all 'f-something' instead of 'something-f'
            new_contrast_str = ['f-' contrasts{cc}(1:2)]; % flip the string around 
            group_data_curr = MRIread([group_dir hemis{hh} '.ces.localizer_groupavg_' new_contrast_str '.glmres/osgm/z.mgh']);
            group_data_curr.vol = -group_data_curr.vol; % flip data to get back to 'something-fixation'
        else
            group_data_curr = MRIread([group_dir hemis{hh} '.ces.localizer_groupavg_' contrasts{cc} '.glmres/osgm/z.mgh']);
        end
        group_data_ROI{cc} = group_data_curr.vol(ROIs{hh}.ind+1); % index numbers in labels/patches start at 0, but indexing starts at 1 in matlab, so add one when indexing like this
    end

    % Lower bound stats at 0 (more interpretable when taking the difference of 2 contrasts)
    group_data_ROI{1}(group_data_ROI{1}<0) = 0;
    group_data_ROI{2}(group_data_ROI{2}<0) = 0;

    % Take difference of 2 contrasts
    group_data_diff{hh} = group_data_ROI{2}' - group_data_ROI{1}';

    % Find axis of greatest change and rotate coordinates so X axis is greatest change axis
    [ROIs{hh}.x, ROIs{hh}.y] = group_find_axis(axis_method, group_data_diff{hh}, ROIs{hh}.x', ROIs{hh}.y', hemi);
end

%% Loop over subjects and fit boundary and gradient models to contrast data
% Lots of stats/diagnostic storage variables
ROI_Ts = cell(N_subjs, N_contrasts, 2); % T stats for each ROI/subj/hemi
winning_model = nan(N_subjs, num_hemis); % which model had best BIC (1 is step, 2 is linear, 3 is hinge)
model_comp = nan(N_subjs, num_hemis); % BIC difference values between winner and next best
winning_rsquare = nan(N_subjs, num_hemis); % r-square of winning model
linear_xdist = nan(N_subjs, num_hemis); % distance of hinge in hinge model (mm)
BICs = nan(N_subjs, num_hemis, 3); % Raw BICs for each model
hinge_direction = nan(N_subjs, num_hemis); % 1 for expected direction (increase along X axis), else 0

% Parameters for how to rotate axis when fitting
deg_step = 1; % degrees
deg_tolerance = 45; % searching +/- 45 degrees off of group axis of greatest change
angles = 0-deg_tolerance:deg_step:0+deg_tolerance;
winning_angle = nan(N_subjs, num_hemis, 3); % which angle was optimal

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

        for cc = 1:length(contrasts)
            % Get contrast t-stat data
            contrast = contrasts{cc};
            tstats = MRIread([subjdir subjCode '/localizer/localizer_contrasts_0sm_' hemi '/' contrast '/t.nii.gz']);
            ROI_tstat = tstats.vol(patch_inds+1); % index numbers in labels/patches start at 0, but indexing starts at 1 in matlab, so add one when indexing like this
            ROI_tstat(ROI_tstat<0) = 0; % Clip negative t-stats at zero so that subtraction of 2 contrasts makes sense
            ROI_Ts{ss,cc,hh} = ROI_tstat;
        end

        % Get difference between contrasts in ROI
        ROI_t_diffs = ROI_Ts{ss,2,hh}  - ROI_Ts{ss,1,hh};
        
        if plot_fits % plot T-stat contrasts in patch before fixing outliers
            figure; subplot(2,1,1);
            scatter(subjROI_x, subjROI_y, [], ROI_t_diffs,'filled');
            clim_set = max(abs(prctile(ROI_t_diffs, [10,90])));
            xlabel('x'); ylabel('y'); cb = colorbar; colormap(redbluedark); clim([-clim_set, clim_set]); ylabel(cb, 'T-stat', 'rotation', 270);
            xlim([min(ROIs{hh}.x),max(ROIs{hh}.x)]); ylim([min(ROIs{hh}.y),max(ROIs{hh}.y)]);
        end

        % Replace outliers with 3 stds away from median value (clipping)
        ROI_t_diffs = filloutliers(ROI_t_diffs, 'clip', 'median'); % 3 stds from median is outlier

        if plot_fits % plot T-stat contrasts in patch after fixing outliers
            subplot(2,1,2);
            scatter(subjROI_x, subjROI_y, [], ROI_t_diffs,'filled');
            clim_set = max(abs(prctile(ROI_t_diffs, [10,90])));
            xlabel('x'); ylabel('y'); cb = colorbar; colormap(redbluedark); clim([-clim_set, clim_set]); ylabel(cb, 'T-stat', 'rotation', 270);
            xlim([min(ROIs{hh}.x),max(ROIs{hh}.x)]); ylim([min(ROIs{hh}.y),max(ROIs{hh}.y)]);
        end
        
        % Rename for clarity/ease
        xs = subjROI_x;
        ys = subjROI_y;
        Ts = ROI_t_diffs';
        group_data_curr = group_data_diff{hh}(patch_inds_mask); % limit group data to inside subj ROI

        % Start plot
        if plot_fits % plot x-coords and T vals
            f1 = figure;
            subplot(2,2,1);
            scatter(xs, Ts);
            title('original');
        end

        %% Fit step-wise function to data (boundary model)
        [gof_step, info_step, winning_angle(ss,hh,1), xs_step, ys_step, fit_res_step] = fit_GB_model('step', xs, ys, Ts, group_data_curr, angles);

        if plot_fits % plot step model fit
            subplot(2,2,2);
            scatter3(xs_step, ys_step, Ts); xlabel('x'); ylabel('y'); zlabel('t-stats'); hold on; %scatter(xs_step, step_results_best);
            plot(fit_res_step)
            view([0 0]);
            title(['step model | rotate ' num2str(winning_angle(ss,hh,1)) ' | rsqr: ' num2str(round(gof_step.rsquare,3))])
        end

        %% Fit 2D plane to data (linear gradient model)
        if strcmp(axis_choice, 'step') % use step model axis 
            angles_linear = winning_angle(ss,hh,1);
        else % find best axis for linear fit
            angles_linear = angles;
        end

        [gof_linear, info_linear, winning_angle(ss,hh,2), xs_linear, ys_linear, fit_res_linear] = ...
            fit_GB_model('linear', xs, ys, Ts, group_data_curr, angles_linear);

        if plot_fits % plot fit for linear model
            subplot(2,2,3);
            scatter3(xs_linear, ys_linear, Ts); xlabel('x'); ylabel('y'); zlabel('t-stats'); hold on; %
            plot(fit_res_linear);
            view([0 0]);
            title(['linear model | rotate ' num2str(winning_angle(ss,hh,2)) ' | rsqr: ' num2str(round(gof_linear.rsquare,3))])
        end

        %% Fit 2D plane to subsection of data (hinge function gradient model)
        if strcmp(axis_choice, 'step') % use step model axis
            angles_hinge = winning_angle(ss,hh,1);
        else % find best axis for hinge fit
            angles_hinge = angles;
        end

        [gof_hinge, info_hinge, winning_angle(ss,hh,3), xs_hinge, ys_hinge, fit_res_hinge] = ...
            fit_GB_model('hinge', xs, ys, Ts, group_data_curr, angles_hinge);

        hinge_direction(ss,hh) = fit_res_hinge.a < fit_res_hinge.b; % see if hinge direction is increase or decrease

        if plot_fits % plot hinge model fit
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

        % Compare BIC between models
        [~, bic] = aicbic([LL_step; LL_linear; LL_hinge], [info_step.numparam; info_linear.numparam; info_hinge.numparam]);
        [~,ind] = min(bic); % get ind of min BIC model
        disp([models{ind}]) % print winning model name
        winning_model(ss,hh) = ind;
        model_comp(ss,hh) = bic(ind) - min(bic(~ismember(1:3,ind)));  % compare best model BIC to next best model BIC
        BICs(ss,hh,:) = bic; % record all BICs

        % Check if the winning model actually fit the data well
        gofs = {gof_step, gof_linear, gof_hinge};
        disp(['winning model r^2: ' num2str(round(gofs{ind}.rsquare,3))]);
        winning_rsquare(ss,hh) = gofs{ind}.rsquare;

        % If hinge model is the winner, check whether the linear piece is large enough to cross multiple voxels
        if ind==3
            linear_xdist(ss,hh) = fit_res_hinge.x2 - fit_res_hinge.x1; % just take difference between end of hinge and beginning of hinge on x axis 
        end

        if plot_fits % make overall title for 2x2 fit plot
            sgtitle([subjCode ' ' hemi ' | Winner:' num2str(winning_model(ss,hh)) ' | BIC Diff:' num2str(round(model_comp(ss,hh),2)) ' | r^2:' num2str(round(winning_rsquare(ss,hh),3)) ' | hingedist: ' num2str(round(linear_xdist(ss,hh),2))])
        end

    end
end

%% Basic counts to compare models
sum(winning_model==3) % how many times did hinge model win?
sum(winning_model==3 & hinge_direction==1) % in the hypothesized direction?
sum(winning_model==3 & winning_rsquare>0.1) % with rsquared over 0.1?
sum(winning_model==3 & model_comp<-10) % more than 10 BIC than next best model?
sum(winning_model==3 & linear_xdist>4) % more than 4mm hinge length?

total_good = sum(winning_model==3 & winning_rsquare>0.1 & model_comp<-10 & linear_xdist>4 & hinge_direction==1) % put them all together

%% Group plots/tests
BIC_allhemis = reshape(BICs, [size(BICs,1)*size(BICs,2), size(BICs,3)]); % collapse hemispheres
hingedir_allhemis = reshape(hinge_direction, [size(hinge_direction,1)*size(hinge_direction,2),1]); % collapse hemispheres
winning_model_allhemis = reshape(winning_model, [size(winning_model,1)*size(winning_model,2),1]); % collapse hemispheres

% Plot BIC differences
new_N = (N_subjs*2) - sum(hingedir_allhemis==0 | isnan(hingedir_allhemis)); % calculate number of subjs/hemis with hinge in correct direction
BIC_diffs = BIC_allhemis(hingedir_allhemis==1,1)-BIC_allhemis(hingedir_allhemis==1,3); % when hinge was in correct direction, what was the BIC difference bt hinge and step?
m = mean(BIC_diffs);
figure; 
subplot(1,3,1);
swarmchart(ones(new_N,1), BIC_diffs, [], 'c', 'filled'); hold on;
swarmchart(1, m, 100, 'r', 'filled');
yline(10, '--r'); ylim([min(BIC_diffs)-20, max(BIC_diffs)+20]);
ylabel('BIC Difference');
xlim([-2,4]); xticks([]);
title({['Step vs Hinge Model | \mu=' num2str(round(m,2))]})

% Plot hinge distances
xdists_hingewinners = linear_xdist(winning_model==3 & hinge_direction==1);
subplot(1,3,2);
swarmchart(ones(length(xdists_hingewinners)), xdists_hingewinners, [], 'c', 'filled'); hold on;
swarmchart(1, mean(xdists_hingewinners), 100, 'r', 'filled');
yline(4, '--r');
xlim([-2,4]); xticks([]);
ylabel('Hinge Length Across Cortex (mm)');
title({['Hinge Length | \mu=' num2str(round(mean(xdists_hingewinners),2))]})

% Plot hinge R squareds
rsqr_hingewinners = winning_rsquare(winning_model==3 & hinge_direction==1);
subplot(1,3,3);
swarmchart(ones(length(rsqr_hingewinners)), rsqr_hingewinners, [], 'c', 'filled'); hold on;
swarmchart(1, mean(rsqr_hingewinners), 100, 'r', 'filled');
yline(0.1, '--r');
xlim([-2,4]); xticks([]);
ylabel('R-squared');
title(['Hinge Model R-squared | \mu=', num2str(round(mean(rsqr_hingewinners),2))])

sgtitle([replace(ROI_name,'_','-') ' ' modailty ' | N = ' num2str(new_N) ', ' num2str(length(rsqr_hingewinners))]);

%% Helper functions %%
function LL= loglikelihood(num_obs, rmse, sse)
    LL = -0.5 * num_obs * log( 2*pi*(rmse^2) ) - (1/(2*(rmse^2))) * sse;
end

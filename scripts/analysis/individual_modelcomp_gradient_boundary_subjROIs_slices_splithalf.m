%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to compare a gradient-based and
%%% boundary-based models of change in function (contrast t-stats) across an ROI in
%%% individual subjects. Split half analysis
%%%
%%% Tom Possidente - Septemeber 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables
hemis = {'lh', 'rh'};
num_hemis = length(hemis);
fs_num = 163842; % number of vertices in fsaverage

% Get participant IDs
subjCodes = {'MK','AB', 'AD','LA','AE','TP','NM','AF','AG','AI','GG','UV','PQ','KQ','LN','RT','PT','PL','NS','SL'};
N_subjs = length(subjCodes);

% Set up directories
subjdir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';
label_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/probabilistic_allROIs/';
group_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/grouplevel/';
subj_ROIs = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/subj_specific_01/';

% Set models/methods used
contrasts = {'vP-f', 'vA-vP'}; % which functional data contrasts to use (supramodal is 'aPvP-f', 'vAaA-vPaP'). P = passive/sensorimotor control, A = active/WM, v = visual, a = auditory, f = fixation
N_contrasts = length(contrasts);
modality = 'visual'; % visual, auditory, supramodal, visaud
ROI_name = 'preSMA'; % Which ROI to look at: preSMA, inf_lat_frontal, aINS, midIFS, sup_lat_frontal, aIPS
models = {'step', 'hinge'};
axis_method = 'average'; % regression (uses regression to find axis of largest difference) or average (uses weighted average of positive and negative pts to make line)
axis_choice = 'step'; % step (use step function axis for all models) or individual (use best axis for each model individually)

plot_fits = false; % plot out individual subj fits (debugging only unless you want a ~100 plots)

dist_thresh = 4.8; % Indicates (in mm) minimum distance ramp section of hinge function must be to be considered a gradient instead of a boundary

%% Load probabilistic ROI (flat patch) to use for determining group level axis of greatest change
ROI_lh = read_patch([label_dir hemis{1} '.' ROI_name '_prob_thresh5_flat.patch']);
ROI_rh = read_patch([label_dir hemis{2} '.' ROI_name '_prob_thresh5_flat.patch']);
ROIs = {ROI_lh, ROI_rh};

%% Compute group level axis of greatest change (used to initialize direction of sensory-WM functional change)
group_data_diff = cell(2,1);
group_angle = nan(2,1);
percent_neg = nan(2,1);

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
    percent_neg(hh) = sum(group_data_ROI{1}<0 | group_data_ROI{2}<0) / length(group_data_ROI{1});

    % Take difference of 2 contrasts
    group_data_diff{hh} = group_data_ROI{2}' - group_data_ROI{1}';

    % Find axis of greatest change and rotate coordinates so X axis is greatest change axis
    [ROIs{hh}.x, ROIs{hh}.y, group_angle(hh)] = group_find_axis(axis_method, group_data_diff{hh}, ROIs{hh}.x', ROIs{hh}.y', hemi);
end

%save(['flatpatch_ROI_rotated_' ROI_name '.mat'], 'ROIs', 'group_angle');

%% Loop over subjects and fit boundary and gradient models to contrast data
% Lots of stats/diagnostic storage variables
Ts_store = cell(N_subjs, N_contrasts, 2, 2); % T stats for each ROI/subj/hemi/oddeven
winning_model = nan(N_subjs, num_hemis); % which model had best BIC (1 is step, 2 is linear, 3 is hinge)
winning_rsquare = nan(N_subjs, num_hemis); % r-square of winning model
linear_xdist = nan(N_subjs, num_hemis); % distance of hinge in hinge model (mm)
change_direction = nan(N_subjs, num_hemis, length(models)); % 1 for expected direction (increase along X axis), else 0
comp_metric = nan(N_subjs, num_hemis, length(models));

% Parameters for how to rotate axis when fitting
deg_step = 1; % degrees
deg_tolerance = 60; % searching +/- deg_toloerance degrees off of group axis of greatest change
angles = 0-deg_tolerance:deg_step:0+deg_tolerance;
perc_neg = nan(N_subjs, num_hemis, N_contrasts);

tic;
for hh = 1:num_hemis
    hemi = hemis{hh};

    for ss = 1:N_subjs
        subjCode = subjCodes{ss};

        % Get subj ROI
        subjROI_path = [subj_ROIs hemi '.' ROI_name '_' modality '_' subjCode '.label'];
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
            % odd runs
            tstats = MRIread([subjdir subjCode '/localizer/localizer_contrasts_0sm_' hemi '_oddruns/' contrast '/t.nii.gz']);
            ROI_tstat = tstats.vol(patch_inds+1); % index numbers in labels/patches start at 0, but indexing starts at 1 in matlab, so add one when indexing like this
            ROI_tstat(ROI_tstat<0) = 0; % Clip negative t-stats at zero so that subtraction of 2 contrasts makes sense
            Ts_store{ss,cc,hh,1} = ROI_tstat;
            %even runs
            tstats = MRIread([subjdir subjCode '/localizer/localizer_contrasts_0sm_' hemi '_evenruns/' contrast '/t.nii.gz']);
            ROI_tstat = tstats.vol(patch_inds+1); % index numbers in labels/patches start at 0, but indexing starts at 1 in matlab, so add one when indexing like this
            ROI_tstat(ROI_tstat<0) = 0; % Clip negative t-stats at zero so that subtraction of 2 contrasts makes sense
            Ts_store{ss,cc,hh,2} = ROI_tstat;
        end

        % Get difference between contrasts in ROI
        ROI_t_diffs = squeeze(cell2mat(Ts_store(ss,2,hh,:)) - cell2mat(Ts_store(ss,1,hh,:)));

        % Replace outliers with 3 stds away from median value (clipping)
        ROI_t_diffs = filloutliers(ROI_t_diffs, 'clip', 'median'); % 3 stds from median is outlier

        % Rename for clarity/ease
        xs = subjROI_x;
        ys = subjROI_y;
        Ts = ROI_t_diffs';
        group_data_curr = group_data_diff{hh}(patch_inds_mask); % limit group data to inside subj ROI

        % Start plot
        if plot_fits % plot x-coords and T vals
            f1 = figure;
            subplot(2,2,1);
            scatter(xs, ys, [], Ts(1,:)', 'filled');
            clim_set = max(abs(prctile(Ts(:), [10,90])));
            xlabel('x'); ylabel('y'); cb = colorbar; colormap(redbluedark); clim([-clim_set, clim_set]); ylabel(cb, 'T-stat', 'rotation', 270);
            xlim([min(xs),max(xs)]); ylim([min(ys),max(ys)]);
            title('original (odd)');
            f2 = figure;
            subplot(2,2,1);
            scatter(xs, ys, [], Ts(2,:)', 'filled');
            clim_set = max(abs(prctile(Ts(:), [10,90])));
            xlabel('x'); ylabel('y'); cb = colorbar; colormap(redbluedark); clim([-clim_set, clim_set]); ylabel(cb, 'T-stat', 'rotation', 270);
            xlim([min(xs),max(xs)]); ylim([min(ys),max(ys)]);
            title('original (even)');
        end

        %% Fit step function to data (boundary model)
        % First find angle that maximizes step model fit for whole ROI (no slices)
        [winning_angle_odd, ~]= individual_find_axis_noslice('step', xs, ys, Ts(1,:), group_data_curr, angles);
        [winning_angle_even, ~]= individual_find_axis_noslice('step', xs, ys, Ts(2,:), group_data_curr, angles);

        % Fit at winning angle per slice
        [metrics_step_odd, ~, ~, xs_new_odd, ys_new_odd, ~, bins_step_odd] = fit_GB_model_slices_splithalf_CV('step', xs, ys, Ts(1,:), group_data_curr, winning_angle_odd, Ts(2,:));
        [metrics_step_even, ~, ~, xs_new_even, ys_new_even, ~, bins_step_even] = fit_GB_model_slices_splithalf_CV('step', xs, ys, Ts(2,:), group_data_curr, winning_angle_even, Ts(1,:));

        if plot_fits % plot step model fit
            figure(f1);
            subplot(2,2,2);
            scatter(xs_new_odd, ys_new_odd, [], Ts(1,:), 'filled');
            clim_set = max(abs(prctile(Ts(:), [10,90])));
            xlabel('x'); ylabel('y'); cb = colorbar; colormap(redbluedark); clim([-clim_set, clim_set]); ylabel(cb, 'T-stat', 'rotation', 270);
            xlim([min(xs_new_odd),max(xs_new_odd)]); ylim([min(ys_new_odd),max(ys_new_odd)]);
            hold on;
            for bb = 1:length(bins_step_odd)-1
                midbin = (bins_step_odd(bb+1) + bins_step_odd(bb)) / 2;
                yline(midbin, ':');
                scatter(metrics_step_odd.p1(bb), midbin, 100, 'g', '^', 'filled')
            end
            title(['step: rotate= ' num2str(winning_angle_odd) ', rsqr=' num2str(round(mean(metrics_step_odd.rsquare, 'omitnan', "Weights", metrics_step_odd.numobs),2))])
            
            figure(f2);
            subplot(2,2,2);
            scatter(xs_new_even, ys_new_even, [], Ts(2,:), 'filled');
            clim_set = max(abs(prctile(Ts(:), [10,90])));
            xlabel('x'); ylabel('y'); cb = colorbar; colormap(redbluedark); clim([-clim_set, clim_set]); ylabel(cb, 'T-stat', 'rotation', 270);
            xlim([min(xs_new_even),max(xs_new_even)]); ylim([min(ys_new_even),max(ys_new_even)]);
            hold on;
            for bb = 1:length(bins_step_even)-1
                midbin = (bins_step_even(bb+1) + bins_step_even(bb)) / 2;
                yline(midbin, ':');
                scatter(metrics_step_even.p1(bb), midbin, 100, 'g', '^', 'filled')
            end
            title(['step: rotate= ' num2str(winning_angle_even) ', rsqr=' num2str(round(mean(metrics_step_even.rsquare, 'omitnan', "Weights", metrics_step_even.numobs),2))])
        end


        %% Fit hinge function gradient model
        [metrics_hinge_odd, ~, ~, xs_new_odd, ys_new_odd, ~] = fit_GB_model_slices_splithalf_CV('hinge', xs, ys, Ts(1,:), group_data_curr, winning_angle_odd, Ts(2,:)); 
        [metrics_hinge_even, ~, ~, xs_new_even, ys_new_even, ~] = fit_GB_model_slices_splithalf_CV('hinge', xs, ys, Ts(2,:), group_data_curr, winning_angle_even, Ts(1,:));

        if plot_fits % plot
            figure(f1);
            subplot(2,2,4);
            scatter(xs_new_odd, ys_new_odd, [], Ts(1,:), 'filled');
            clim_set = max(abs(prctile(Ts(:), [10,90])));
            xlabel('x'); ylabel('y'); cb = colorbar; colormap(redbluedark); clim([-clim_set, clim_set]); ylabel(cb, 'T-stat', 'rotation', 270);
            xlim([min(xs_new_odd),max(xs_new_odd)]); ylim([min(ys_new_odd),max(ys_new_odd)]);
            hold on;
            for bb = 1:length(bins_step_odd)-1
                midbin = (bins_step_odd(bb+1) + bins_step_odd(bb)) / 2;
                yline(midbin, ':');
                scatter(metrics_hinge_odd.p1(bb), midbin, 100, 'g', '^', 'filled')
                scatter(metrics_hinge_odd.p1(bb)+metrics_hinge_odd.p2(bb), midbin, 100, 'g', '>', 'filled')
            end
            title(['hinge: rsqr=' num2str(round(mean(metrics_hinge_odd.rsquare, 'omitnan', "Weights", metrics_step_odd.numobs),2))])  

            figure(f2);
            subplot(2,2,4);
            scatter(xs_new_even, ys_new_even, [], Ts(2,:), 'filled');
            clim_set = max(abs(prctile(Ts(:), [10,90])));
            xlabel('x'); ylabel('y'); cb = colorbar; colormap(redbluedark); clim([-clim_set, clim_set]); ylabel(cb, 'T-stat', 'rotation', 270);
            xlim([min(xs_new_even),max(xs_new_even)]); ylim([min(ys_new_even),max(ys_new_even)]);
            hold on;
            for bb = 1:length(bins_step_even)-1
                midbin = (bins_step_even(bb+1) + bins_step_even(bb)) / 2;
                yline(midbin, ':');
                scatter(metrics_hinge_even.p1(bb), midbin, 100, 'g', '^', 'filled')
                scatter(metrics_hinge_even.p1(bb)+metrics_hinge_even.p2(bb), midbin, 100, 'g', '>', 'filled')
            end
            title(['hinge: rsqr=' num2str(round(mean(metrics_hinge_even.rsquare, 'omitnan', "Weights", metrics_step_even.numobs),2))])    
        end

        %% Take weighted means of metrics
        step_means_odd = mean(metrics_step_odd, 'omitnan', 'Weights', metrics_step_odd.numobs);
        step_means_even = mean(metrics_step_even, 'omitnan', 'Weights', metrics_step_even.numobs);

        hinge_means_odd = mean(metrics_hinge_odd, 'omitnan', 'Weights', metrics_step_odd.numobs);
        hinge_means_even = mean(metrics_hinge_even, 'omitnan', 'Weights', metrics_step_even.numobs);

        %% Take means of odd and even splits
        step_means = mean([step_means_odd; step_means_even]);
        hinge_means = mean([hinge_means_odd; hinge_means_even]);

        %% Get hinge distance
        linear_xdist(ss,hh) = hinge_means.p2;

        %% Store CV correlation comparison metric
        comp_metric(ss,hh,:) = [step_means.CV_corr, hinge_means.CV_corr];

        %% Get winning model (highest corr)
        [~, winning_model(ss,hh)] = max([step_means.CV_corr, hinge_means.CV_corr]);

        %% Check if the winning model actually fit the data well
        rsquares = [mean([step_means_odd.rsquare, step_means_even.rsquare]), mean([hinge_means_odd.rsquare, hinge_means_even.rsquare])];
        disp(['winning model r^2: ' num2str(round(rsquares(winning_model(ss,hh)),3))]);
        winning_rsquare(ss,hh) = rsquares(winning_model(ss,hh));

        %% Check direction of change for each model
        change_direction(ss,hh,1) = round(mean(double(step_means.p2 < step_means.p3), 'omitnan', 'Weights', step_means.numobs));
        change_direction(ss,hh,2) = round(mean(double(hinge_means.p3 < hinge_means.p4), 'omitnan', 'Weights', hinge_means.numobs));        

    end
end
toc

%% Check the change direction of each winning model
good_winning_direction = zeros(N_subjs,num_hemis);
for mm = 1:length(models)
    good_winning_direction = good_winning_direction + ((winning_model==mm) .* change_direction(:,:,mm));
end

%% CV model comparison
% Reject bad winning directions
comp_metric(change_direction==0) = nan;
linear_xdist(change_direction(:,:,2)==0) = nan;

% take mean across hemispheres
comp_metric = squeeze(mean(comp_metric, 2, 'omitnan'));

% Get means and SEs 
mean_step_CVcorr = mean(comp_metric(:,1),'omitnan');
SE_step_CVcorr = std(comp_metric(:,1),'omitnan')/sqrt(sum(~isnan(comp_metric(:,1))));
mean_hinge_CVcorr = mean(comp_metric(:,2),'omitnan');
SE_hinge_CVcorr = std(comp_metric(:,2),'omitnan')/sqrt(sum(~isnan(comp_metric(:,2))));

disp(['Mean boundary CV corr: ' num2str(round(mean_step_CVcorr,3))])
disp(['SE boundary CV corr: ' num2str(round(SE_step_CVcorr,3))])
disp(['Mean gradient CV corr: ' num2str(round(mean_hinge_CVcorr,3))])
disp(['SE gradient CV corr: ' num2str(round(SE_hinge_CVcorr,3))])

% Correlation comparison
[h,p,ci,stats] = ttest(comp_metric(:,1), comp_metric(:,2)) % paired t-test

save(['CV_modelcomp_results_' ROI_name '_' modality '_nobaddir.mat'], 'mean_step_CVcorr', 'SE_step_CVcorr', 'mean_hinge_CVcorr', 'SE_hinge_CVcorr', 'h','p','ci','stats', 'linear_xdist');




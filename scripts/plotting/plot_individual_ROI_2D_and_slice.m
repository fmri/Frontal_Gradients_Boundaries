%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to plot one
%%% subj/modality/hemisphere/ROI's 2D flatmap with contrast visualization.
%%%
%%% Tom Possidente - September 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables

modality = 'visual';
ROI_name = 'sup_lat_frontal';
subjCodes = {'AB'};  % subjCodes = {'MK','AB', 'AD','LA','AE','TP','NM','AF','AG','AI','GG',
                               % 'UV','PQ','KQ','LN','RT','PT','PL','NS','SL'};
N_subjs = length(subjCodes);
hemi = 'lh';
slice_ind = [11];

save_png = true;

subj_ROIs = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/subj_specific_01/';
label_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/probabilistic_allROIs/';
subjdir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';
angle_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/Individual_ROI_models/';

if strcmp(modality, 'visual')
    contrasts = {'vP-f', 'vA-vP'};
elseif strcmp(modality, 'auditory')
    contrasts = {'aP-f', 'aA-aP'};
else
    error('unrecognized modality')
end

if strcmp(hemi, 'lh')
    hh = 1;
elseif strcmp(hemi, 'rh')
    hh = 2;
else
    error('unrecognized hemisphere')
end


ROI_lh = read_patch([label_dir 'lh.' ROI_name '_prob_thresh5_flat.patch']);
ROI_rh = read_patch([label_dir 'rh.' ROI_name '_prob_thresh5_flat.patch']);
ROIs = {ROI_lh, ROI_rh};

Ts_store = cell(2,1);

slice_size = 2.2; %mm

for nn = 1:N_subjs
    %% Get data
    subjROI_path = [subj_ROIs hemi '.' ROI_name '_' modality '_' subjCodes{nn} '.label'];
    if ~isfile(subjROI_path) % if this ROI doesn't exist, the subj didn't have enough good signal/vertices in this ROI, skip it
        error('ROI not found');
    end
    
    subjROI = readtable(subjROI_path, 'FileType','text'); % load label file for individual subj ROI
    patch_inds_mask = ismember(ROIs{hh}.ind, subjROI{:,1}); % Get flat patch indices/x/y for each vertex in subj ROI
    patch_inds = ROIs{hh}.ind(patch_inds_mask);
    subjROI_x = ROIs{hh}.x(patch_inds_mask);
    subjROI_y = ROIs{hh}.y(patch_inds_mask);
    
    for cc = 1:length(contrasts)
        % Get contrast t-stat data
        contrast = contrasts{cc};
        tstats = MRIread([subjdir subjCodes{nn} '/localizer/localizer_contrasts_0sm_' hemi '/' contrast '/t.nii.gz']);
        ROI_tstat = tstats.vol(patch_inds+1); % index numbers in labels/patches start at 0, but indexing starts at 1 in matlab, so add one when indexing like this
        ROI_tstat(ROI_tstat<0) = 0; % Clip negative t-stats at zero so that subtraction of 2 contrasts makes sense
        Ts_store{cc} = ROI_tstat;
    end
    
    % Get difference between contrasts in ROI
    ROI_t_diffs = Ts_store{2}  - Ts_store{1};
    
    % Replace outliers with 3 stds away from median value (clipping)
    ROI_t_diffs = filloutliers(ROI_t_diffs, 'clip', 'median'); % 3 stds from median is outlier
    
    % Rename for clarity/ease
    xs = subjROI_x;
    ys = subjROI_y;
    Ts = ROI_t_diffs';
    
    % Rotate to winning angle (SOMETHING WRONG HERE? NOT ROTATING AS EXPECTED)
    angle = load([angle_dir ROI_name '_' modality '_winning_angles.mat']);
    theta = deg2rad(angle.winning_angles(ismember(angle.subjCodes, subjCodes{nn}), ismember(angle.hemis, hemi)));
    xs_new = xs * cos(theta) - ys * sin(theta); % get new rotated x and y coords
    ys_new = xs * sin(theta) + ys * cos(theta);
    
    % plot
    f1 = figure;
    scatter(xs_new, ys_new, [], Ts, 'filled');
    clim_set = max(abs(prctile(Ts, [10,90])));
    xlabel('x'); ylabel('y'); cb = colorbar; colormap(redbluedark); clim([-clim_set, clim_set]); ylabel(cb, 'T-stat', 'rotation', 270);
    middle_x = mean([max(xs_new), min(xs_new)]); 
    middle_y = mean([max(ys_new), min(ys_new)]); 
    xlim([middle_x-20, middle_x+20]); ylim([middle_y-20, middle_y+20]);
    title([subjCodes{nn} ' ' modality ' ' hemi ' ' replace(ROI_name,'_', ' ')]);
    fontsize(16, 'points');
    
    bins = min(ys_new):slice_size:max(ys_new);
    n_slices = length(bins) - 1;
    
    for bb = 1:n_slices
        midbin = (bins(bb+1) + bins(bb)) / 2;
        yline(midbin, ':');
    end
    f1.Position = [496   279   694   500];
    
    if save_png
        print('-dpng', '-r600', ['2Dcontrast_flatmap_' subjCodes{nn} '_' modality '_' ROI_name '_' hemi '.png'])
    end
    
    
    %% get models for specific slice
    slice_bin = bins(slice_ind(nn):slice_ind(nn)+1);
    slice_mask = ys_new>=slice_bin(1) & ys_new<slice_bin(2);
    slice_xs = xs_new(slice_mask);
    slice_ys = ys_new(slice_mask);
    slice_Ts = Ts(slice_mask);
    
    slice_group_data = [];
    
    if length(slice_Ts) < 20
        error(['Fewer than 20 data points in slice ' num2str(slice_ind(nn)) '/' num2str(n_slices)])
    elseif all(slice_Ts==0)
        error(['All data points equal zero in slice ' num2str(slice_ind(nn)) '/' num2str(n_slices)])
    end
    
    step_model = fittype('x2*(x<x1) + x3*(x>=x1)',... % step function model
        'dependent', 'z',...
        'independent', {'x'}, ...
        'coefficients', {'x1','x2', 'x3'}); % x1 is step location, x2 is pre-step, x3 is post-step
    
    startpoint_guesses = [mean(slice_xs), -1.5, 1.5]; % if no group data, hardcoded guesses
    [fit_res_step, gof_step, info_step] = fit(slice_xs', slice_Ts, step_model, 'StartPoint', startpoint_guesses); % actually fit
    
    model_hinge = fittype('a*(x < x1) + b*(x > (x1+x2)) + ( a + ( (b-a)/((x1+x2)-x1)  * (x-x1) ) ) * ( (x > x1) & (x < (x1+x2)) )', ... % piecewise function, constant from -Inf to x1, linear from x1 to x2, constant from x2 to Inf. Must include y in the function even if it has no mathematical effect
        'dependent', 'z',...
        'independent', {'x'}, ...
        'coefficients', {'a','b','x1', 'x2'}); % a is z from -Inf to x1, b is z from x2 to Inf, x1 is where to start linear piece, x2 is where to end linear piece
    
    startpoint_guesses = [-1.5, 1.5, prctile(slice_xs, 25), prctile(slice_xs, 50)]; 
    
    lower_bounds = [-10, -10, min(slice_xs), 0];
    upper_bounds = [10, 10, max(slice_xs), max(slice_xs)-min(slice_xs)];
    
    [fit_res_hinge, gof_hinge, info_hinge] = fit(slice_xs', slice_Ts, model_hinge, 'StartPoint', startpoint_guesses,...
        'Lower', lower_bounds,...
        'Upper', upper_bounds); 
    
    
    % Plot slice fits
    f2 = figure;
    subplot(1,2,1);
    scatter(slice_xs', slice_Ts, 'filled'); hold on;
    plot([min(slice_xs), fit_res_step.x1, fit_res_step.x1, max(slice_xs)], ...
        [fit_res_step.x2, fit_res_step.x2, fit_res_step.x3, fit_res_step.x3], 'LineWidth',3, 'Color', [0.8 0.7 0])
    title(['r2=' num2str(round(gof_step.rsquare,3))])
    xlim([middle_x-15, middle_x+15]);

    subplot(1,2,2);
    scatter(slice_xs', slice_Ts, 'filled'); hold on;
    plot([min(slice_xs), fit_res_hinge.x1, fit_res_hinge.x1+fit_res_hinge.x2, max(slice_xs)],...
        [fit_res_hinge.a, fit_res_hinge.a, fit_res_hinge.b, fit_res_hinge.b], 'LineWidth',3, 'Color', [0.8 0.7 0]);
    title(['r2=' num2str(round(gof_hinge.rsquare,3)) ' | hinge=' num2str(round(fit_res_hinge.x2,2)) 'mm'])
    xlim([middle_x-15, middle_x+15]);
    fontsize(16, 'points');
    f2.Position = [    1           1        1154         450];
    
    if save_png
        print('-dpng', '-r600', ['1Dcontrast_slice_' num2str(slice_ind(nn)) '_' subjCodes{nn} '_' modality '_' ROI_name '_' hemi '.png'])
    end
end

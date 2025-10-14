%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to compare a gradient-based and
%%% boundary-based model of change in function (contrast PSC) across an ROI in
%%% individual subjects
%%%
%%% TODOs: 
%%% - Technically we should be calculating a separate group average
%%% axis of largest change and initialial fit params using a leave one out
%%% method for each subject 
%%%
%%% Tom Possidente - Septemeber 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/');
ccc;

%% Initialize Key Variables
contrasts = {'aPvP-f', 'vAaA-vPaP'};
latmed = 'lateral_middle'; % lateral or medial

hemis = {'lh', 'rh'};
fs_num = 163842;

experiment_name = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode(~ismember(subjDf_cut.subjCode, {'AH', 'SL', 'RR', 'PP', 'MM'}));
N_subjs = length(subjCodes);

subjdir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';
fs_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/fsaverage/surf/';
label_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/';

models = {'step', 'linear', 'hinge'};
plot_fits = true;

axis_method = 'average'; % regression (uses regression to find axis of largest difference) or average (uses weighted average of positive and negative pts to make line)

%% Load label file for probabilistic ROIs
label_lh = readtable([label_dir hemis{1} '.' latmed '_VisAudWM_combined_TFCE.label'], 'FileType','text');
label_rh = readtable([label_dir hemis{2} '.' latmed '_VisAudWM_combined_TFCE.label'], 'FileType','text');
labels = {label_lh, label_rh};

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
            new_x{hh} = v' / norm(v);  % normalize
            
            % Perpendicular vector (90° rotation)
            new_y{hh} = [-new_x{hh}(2); new_x{hh}(1)];

            % Plot 3D data surface
            figure;
            scatter3(ROIs{hh}.x, ROIs{hh}.y, group_data_diff{hh});
            hold on;
            plot([weighted_pos(1), weighted_neg(1)], [weighted_pos(2), weighted_neg(2)], 'LineWidth',5, 'Color', 'r');
            xlabel('x'); ylabel('y'); zlabel('PSC')
    end
end

%% Rotate lh and rh flatmap coordinates to new x and y axes
for hh = 1:length(hemis)
    xy = [ROIs{hh}.x; ROIs{hh}.y];
    R = [new_x{hh}, new_y{hh}]';
    xy_rotated = R * xy;
    ROIs{hh}.x = xy_rotated(1,:);
    ROIs{hh}.y = xy_rotated(2,:);
    
    % Plot psc data along new axes
    % [xq,yq] = meshgrid(linspace(min(ROIs{hh}.x), max(ROIs{hh}.x), 100),...
    %                linspace(min(ROIs{hh}.y), max(ROIs{hh}.y), 100));
    % zq = griddata(ROIs{hh}.x, ROIs{hh}.y, group_data_diff{hh}, xq, yq, 'linear'); % interpolate to 100x100 grid
    figure;
    % surf(xq, yq, zq); % plot surface    
    scatter3(ROIs{hh}.x, ROIs{hh}.y, group_data_diff{hh});
    xlabel('x'); ylabel('y'); zlabel('PSC');
end


%% Loop over subjects and fit boundary and gradient models to contrast data
ROI_pscs = {nan(length(ROIs{1}.ind),2,N_subjs), nan(length(ROIs{2}.ind),2,N_subjs)};
winning_model = nan(N_subjs, length(hemis));
model_comp = nan(N_subjs, length(hemis));
winning_rsquare = nan(N_subjs, length(hemis));
linear_xdist = nan(N_subjs, length(hemis));
BICs = nan(N_subjs, length(hemis), 3);
hinge_direction = nan(N_subjs, length(hemis));

for hh = 1:length(hemis)
    hemi = hemis{hh};
    
    % Calculate initial guesses for boundary model parameters using group-level data
    x_midpoint = mean(ROIs{hh}.x); % boundary step function location guess
    min_guess = mean(group_data_diff{hh}(group_data_diff{hh}<0)); % boundary model lower step guess
    max_guess = mean(group_data_diff{hh}(group_data_diff{hh}>0)); % boundary model upper step guess
    lower_bound = min(ROIs{hh}.x); % Lower bound on step location
    upper_bound = max(ROIs{hh}.x); % upper bound on step location

    % Calculate initial guesses for gradient model parameters using group-level data
    X = [ROIs{hh}.x; ones(1,length(ROIs{hh}.x))]';
    coefs = X \ group_data_diff{hh}'; % linear regression on only x coordinates
    slope_guess = coefs(1); % 1st param is slope
    intercept_guess = coefs(2); % 2nd is intercept

    % Calculate initial guesses for linear limited model parameters 
    start_line_guess = prctile(ROIs{hh}.x, 25); % lets say it starts at about 25% of x coordinates 
    end_line_guess = prctile(ROIs{hh}.x, 75); % lets say it ends at about 75% of x coordinates

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
        
        % Fit step-wise funciton to data (boundary model)
        step_model = @(params, x) params(2)*(x < params(1)) + params(3)*(x >= params(1)); % param(1) is step location, 2 is pre-step, 3 is post-step
        lossFunction_step = @(params) sum((ROI_psc_diffs' - step_model(params, ROIs{hh}.x)).^2 );
        optParams = fminsearch(lossFunction_step, [x_midpoint, max_guess, min_guess]);
        step_fit = step_model(optParams, ROIs{hh}.x);
        residuals_step = ROI_psc_diffs' - step_fit;
        info_step.numobs = length(ROIs{hh}.x);
        gof_step.rmse = sqrt(mean(residuals_step.^2));
        gof_step.sse = sum(residuals_step.^2);
        gof_step.rsquare = 1 - (gof_step.sse / sum( (ROI_psc_diffs'-mean(ROI_psc_diffs)).^2 ) );
        info_step.numparam = 3;

        % step_model = fittype('a*(x < x0) + b*(x >= x0) + y*0', ... % step function. Must include y in the function even if it has no mathematical effect
        %                     'dependent', 'z',...
        %                     'independent', {'x','y'}, ... % y is not actually used in the function but include it so that we can use all 3D data to fit (not just the x and z coordinates)
        %                     'coefficients', {'a','b','x0'}); % a is pre-step, b is post-step, x0 is step location on x axis
        % [fit_res_step, gof_step, info_step] = fit([ROIs{hh}.x', ROIs{hh}.y'], ROI_psc_diffs, step_model, ...
        %                                          'StartPoint', [max_guess, min_guess, x_midpoint],...
        %                                          'Lower', [-Inf, -Inf, lower_bound],... % 
        %                                          'upper', [Inf, Inf, upper_bound]);  % 
        % disp(['x0=' num2str(fit_res_step.x0)])

        % Fit 2D plane to data (gradient model)
        linear_model = fittype('x*a + b + y*0',... % linear function
                               'dependent', 'z',...
                               'independent', {'x','y'}, ... % y is not actually used in the function but include it so that we can use all 3D data to fit (not just the x and z coordinates)
                               'coefficients', {'a','b'}); % a is slope, b is intercept

        [fit_res_linear, gof_linear, info_linear] = fit([ROIs{hh}.x', ROIs{hh}.y'], ROI_psc_diffs, linear_model, ...
                                                        'StartPoint', [slope_guess, intercept_guess]);  


        % Fit 2D plane to subsection of data (gradient adjusted model)
        hinge_model = fittype('a*(x < x1) + b*(x > x2) + ( a + ( (b-a)/(x2-x1)  * (x-x1) ) ) * ( (x > x1) & (x < x2) ) + y*0', ... % piecewise function, constant from -Inf to x1, linear from x1 to x2, constant from x2 to Inf. Must include y in the function even if it has no mathematical effect
                            'dependent', 'z',...
                            'independent', {'x','y'}, ... % y is not actually used in the function but include it so that we can use all 3D data to fit (not just the x and z coordinates)
                            'coefficients', {'a','b','x1', 'x2'}); % a is z from -Inf to x1, b is z from x2 to Inf, x1 is where to start linear piece, x2 is where to end linear piece

        [fit_res_hinge, gof_hinge, info_hinge] = fit([ROIs{hh}.x', ROIs{hh}.y'], ROI_psc_diffs, hinge_model, ...
                                                                 'StartPoint', [max_guess, min_guess, start_line_guess, end_line_guess],...
                                                                 'Lower', [-Inf, -Inf, lower_bound, lower_bound],...
                                                                 'Upper', [Inf, Inf, upper_bound, upper_bound]);  

        % hinge_model = @(params, x) params(1)*(x < params(3)) + params(2)*(x > params(4)) + ( params(1) + ( (params(2)-params(1))/(params(4)-params(3))  * (x-params(3)) ) ) .* ( (x > params(3)) & (x < params(4)) ); % param(1) is pre-line z val, 2 post-line z val, 3 is line start location, 4 is line end location
        % lossFunction_hinge = @(params) sum((ROI_psc_diffs' - hinge_model(params, ROIs{hh}.x)).^2 );
        % optParams = fminsearch(lossFunction_hinge, [max_guess, min_guess, start_line_guess, end_line_guess]);
        % hinge_fit = hinge_model(optParams, ROIs{hh}.x);
        % residuals_hinge = ROI_psc_diffs' - hinge_fit;
        % info_hinge.numobs = length(ROIs{hh}.x);
        % gof_hinge.rmse = sqrt(mean(residuals_hinge.^2));
        % gof_hinge.sse = sum(residuals_hinge.^2);
        % gof_hinge.rsquare = 1 - (gof_hinge.sse / sum( (ROI_psc_diffs'-mean(ROI_psc_diffs)).^2 ) );
        info_hinge.numparam = 4;
        hinge_direction(ss,hh) = fit_res_hinge.a < fit_res_hinge.b; 

        % Calculate log likelihood of each model
        LL_step = -0.5 * info_step.numobs * log( 2*pi*(gof_step.rmse^2) ) - (1/(2*(gof_step.rmse^2))) * gof_step.sse;
        LL_linear = -0.5 * info_linear.numobs * log( 2*pi*(gof_linear.rmse^2) ) - (1/(2*(gof_linear.rmse^2))) * gof_linear.sse;
        LL_hinge = -0.5 * info_hinge.numobs * log( 2*pi*(gof_hinge.rmse^2) ) - (1/(2*(gof_hinge.rmse^2))) * gof_hinge.sse;

        % Compare AIC/BIC between models
        [aic, bic] = aicbic([LL_step; LL_linear; LL_hinge], [info_step.numparam; info_linear.numparam; info_hinge.numparam]);
        [~,ind] = min(bic);
        disp([models{ind}])
        winning_model(ss,hh) = ind;
        model_comp(ss,hh) = bic(ind) - min(bic(~ismember(1:3,ind)));  % compare best model BIC to next best model BIC
        BICs(ss,hh,:) = bic;

        % Check if the "winning" model fit the data well
        gofs = {gof_step, gof_linear, gof_hinge};
        disp(['winning model r^2: ' num2str(round(gofs{ind}.rsquare,3))]);
        winning_rsquare(ss,hh) = gofs{ind}.rsquare;

        % If linear limited is the winner, check whether the linear piece is large enough to cross multiple voxels
        if ind==3
            mean_y1 = mean(ROIs{hh}.y( (ROIs{hh}.x<fit_res_hinge.x1+0.5) & (ROIs{hh}.x>fit_res_hinge.x1-0.5) ) ); % get mean y coord near x1
            mean_y2 = mean(ROIs{hh}.y( (ROIs{hh}.x<fit_res_hinge.x2+0.5) & (ROIs{hh}.x>fit_res_hinge.x2-0.5) ) ); % get mean y coord near x1
            linear_extent = [fit_res_hinge.x1, mean_y1, ; fit_res_hinge.x2, mean_y2];
            
            % Find vertex nearest (x1, mean_y) and (x2, mean_y) and get RAS coords
            [dist1, vert1_ind] = min(cell2mat(arrayfun(@(x) pdist([linear_extent(1,:); [ROIs{hh}.x(x), ROIs{hh}.y(x)]]), 1:length(ROIs{hh}.x), 'UniformOutput', false) ) );
            if dist1>0.5
                disp('Poor vertex match');
                keyboard;
            end
            vert1 = ROIs{hh}.ind(vert1_ind);
            RAS_vert1 = labels{hh}{labels{hh}.Var1==vert1, 2:4};
            [dist2, vert2_ind] = min(cell2mat(arrayfun(@(x) pdist([linear_extent(2,:); [ROIs{hh}.x(x), ROIs{hh}.y(x)]]), 1:length(ROIs{hh}.x), 'UniformOutput', false) ) );
            if dist2>0.5
                disp('Poor vertex match');
                keyboard;
            end            
            vert2 = ROIs{hh}.ind(vert2_ind);
            RAS_vert2 = labels{hh}{labels{hh}.Var1==vert2, 2:4};

            % Find distance between verts
            linear_xdist(ss,hh) = pdist([RAS_vert1; RAS_vert2]);
        end

                % Plot fit (debugging/visualization)
        if plot_fits
            figure;
            scatter3(ROIs{hh}.x, ROIs{hh}.y, ROI_psc_diffs); xlabel('x'); ylabel('y'); zlabel('psc'); hold on;
            scatter3(ROIs{hh}.x, ROIs{hh}.y, step_fit); hold on;
            %scatter3(ROIs{hh}.x, ROIs{hh}.y, linlim_fit); hold on;
            plot(fit_res_linear);
            plot(fit_res_hinge)
            view([0 0]);
            title([subjCode ' ' hemi ' | Winner:' num2str(winning_model(ss,hh)) ' | BIC Diff:' num2str(model_comp(ss,hh)) ' r^2:' num2str(winning_rsquare(ss,hh)) ])
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


total_good = sum(winning_model==3 & winning_rsquare>0.05 & linear_xdist>4 & hinge_direction==1);
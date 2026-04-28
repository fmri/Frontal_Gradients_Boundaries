function [key_metrics, gofs, model_info, xs_new, ys_new, fit_results, bins] = fit_GB_model_slices(type, xs, ys, Ts, group_data, angle, slice_size)
% FIT_GB_MODEL fits boundary (step function) or gradient (linear, hinge
% functions) to the given data at the supplied angle rotations.
%
%   Note that for all custom functions used with Matlab's "fit" function,
%   the default cost function is non-linear least squares and the default
%   algorithm is trust-region

arguments (Input)
    type {mustBeText, mustBeNonempty} % which model to fit (step, linear, hinge)
    xs (:,1) {mustBeMatrix, mustBeNonNan, mustBeNonempty} % x coordinates
    ys (:,1) {mustBeMatrix, mustBeNonNan, mustBeNonempty} % y coordinates
    Ts (:,1) {mustBeMatrix, mustBeNonNan, mustBeNonempty} % individual data values corresponding to x and y coordinates
    group_data (:,1) {mustBeMatrix, mustBeNonNan, mustBeNonempty} % group data values corresponding to x and y coordinates (used for guessing starting fit values)
    angle (1,1) {mustBeNumeric, mustBeNonNan, mustBeNonempty} % degrees to rotate the coordinates and fit
    slice_size (1,1) {mustBeNumeric} = 4.4 % mm for length of slices in along y axis
end

% Rotate data to desired angle outslide slice loop
theta = deg2rad(angle);
xs_new = xs * cos(theta) - ys * sin(theta); % get new rotated x and y coords
ys_new = xs * sin(theta) + ys * cos(theta);

% Get slice bins
bins = min(ys_new):slice_size:max(ys_new);
n_slices = length(bins) - 1;

% Variable storage initialization
gofs = cell(n_slices,1);
model_info = cell(n_slices,1);
fit_results = cell(n_slices,1);
key_metrics = nan(n_slices, 8); % rsquared, rmse, sse, num_observations, parameters1-4 (some may be blank if model does not have 4 params)

for ss = 1:n_slices

    % Get slice data
    slice_bin = bins(ss:ss+1);
    slice_mask = ys_new>=slice_bin(1) & ys_new<slice_bin(2);
    slice_xs = xs_new(slice_mask);
    slice_ys = ys_new(slice_mask);
    slice_Ts = Ts(slice_mask);
    slice_group_data = group_data(slice_mask);
    if length(slice_Ts) < 20
        disp(['Fewer than 20 data points in slice ' num2str(ss) '/' num2str(n_slices) '... skipping slice']);
        continue
    elseif all(slice_Ts==0)
        disp(['All data points equal zero in slice ' num2str(ss) '/' num2str(n_slices) '... skipping slice']);
        continue
    end

    if strcmp(type, 'step')
        model = fittype('x2*(x<x1) + x3*(x>=x1)',... % step function model
            'dependent', 'z',...
            'independent', {'x'}, ... % y is not actually used in the function but include it so that we can use all 3D data to fit (not just the x and z coordinates)
            'coefficients', {'x1','x2', 'x3'}); % x1 is step location, x2 is pre-step, x3 is post-step

        % Calculate initial guesses for boundary model parameters using group-level data
        startpoint_guesses = [mean(slice_xs), mean([slice_group_data(slice_group_data<0);0]), mean([slice_group_data(slice_group_data>0); 0])]; % estimates for step location, pre-step z value, and post-step z value

        [fit_res_curr, gof_curr, info_curr] = fit(slice_xs, slice_Ts, model, 'StartPoint', startpoint_guesses); % actually fit
        parameters = [fit_res_curr.x1, fit_res_curr.x2, fit_res_curr.x3, nan];
    elseif strcmp(type, 'linear')
        model = fittype('x*a + b',... % linear function
            'dependent', 'z',...
            'independent', {'x'}, ... % y is not actually used in the function but include it so that we can use all 3D data to fit (not just the x and z coordinates)
            'coefficients', {'a','b'}); % a is slope, b is intercept

        % Calculate initial guesses for gradient model parameters using group-level data
        X = [slice_xs, ones(length(slice_xs),1)];
        coefs = X \ slice_group_data; % linear regression on only x coordinates
        slope_guess = coefs(1); % 1st param is slope
        intercept_guess = coefs(2); % 2nd is intercept
        startpoint_guesses = [slope_guess, intercept_guess];

        [fit_res_curr, gof_curr, info_curr] = fit(slice_xs, slice_Ts, model, 'StartPoint', startpoint_guesses); % actually fit
        parameters = [fit_res_curr.a, fit_res_curr.b, nan, nan];
    elseif strcmp(type, 'hinge')
        model = fittype('a*(x < x1) + b*(x > (x1+x2)) + ( a + ( (b-a)/((x1+x2)-x1)  * (x-x1) ) ) * ( (x > x1) & (x < (x1+x2)) )', ... % piecewise function, constant from -Inf to x1, linear from x1 to x2, constant from x2 to Inf. Must include y in the function even if it has no mathematical effect
            'dependent', 'z',...
            'independent', {'x'}, ... % y is not actually used in the function but include it so that we can use all 3D data to fit (not just the x and z coordinates)
            'coefficients', {'a','b','x1', 'x2'}); % a is z from -Inf to x1, b is z from x2 to Inf, x1 is where to start linear piece, x2 is where to end linear piece

        % Calculate initial guesses for gradient model parameters using group-level data
        startpoint_guesses = [max([mean([slice_group_data(slice_group_data<0);0]), min(slice_group_data)]), min([mean([slice_group_data(slice_group_data>0);0]), ...
                              max(slice_group_data)]), prctile(slice_xs, 25), prctile(slice_xs, 50)]; % estimates
        lower_bounds = [-10, -10, min(slice_xs), 0];
        upper_bounds = [10, 10, max(slice_xs), max(slice_xs)-min(slice_xs)];

        [fit_res_curr, gof_curr, info_curr] = fit(slice_xs, slice_Ts, model, 'StartPoint', startpoint_guesses,...
            'Lower', lower_bounds,...
            'Upper', upper_bounds); % actually fit
        parameters = [fit_res_curr.x1, fit_res_curr.x2, fit_res_curr.a, fit_res_curr.b]; 

    else
        error('model type not recognized')
    end

    fit_results{ss} = fit_res_curr;
    gofs{ss} = gof_curr;
    model_info{ss}= info_curr;
    key_metrics(ss,:) = [gof_curr.rsquare, gof_curr.rmse, gof_curr.sse, info_curr.numobs, parameters];
end
key_metrics = array2table(key_metrics, 'VariableNames',{'rsquare', 'rmse', 'sse', 'numobs', 'p1', 'p2', 'p3', 'p4'});
end
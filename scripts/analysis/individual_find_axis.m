function [winning_angle] = individual_find_axis(type, xs, ys, Ts, group_data, angles, slice_size)
% INDIVIDUAL_FIND_AXIS fits boundary (step function) or gradient (linear, hinge
% functions) to the given data at the supplied angle rotations to see which
% is the best angle to use for a given model

arguments (Input)
    type {mustBeText, mustBeNonempty} % which model to fit (step, linear, hinge)
    xs (:,1) {mustBeMatrix, mustBeNonNan, mustBeNonempty} % x coordinates
    ys (:,1) {mustBeMatrix, mustBeNonNan, mustBeNonempty} % y coordinates
    Ts (:,1) {mustBeMatrix, mustBeNonNan, mustBeNonempty} % individual data values corresponding to x and y coordinates
    group_data (:,1) {mustBeMatrix, mustBeNonNan, mustBeNonempty} % group data values corresponding to x and y coordinates (used for guessing starting fit values)
    angles {mustBeMatrix, mustBeNonNan, mustBeNonempty} % degrees to rotate the coordinates and fit (this function will output the info for the best fit out of all angles)
    slice_size (1,1) {mustBeNumeric} = 4.4 % mm for length of slices in along y axis
end

best_rsquare = -100; % initialize

if strcmp(type, 'step')
    model = fittype('x2*(x<x1) + x3*(x>=x1)',... % step function model
        'dependent', 'z',...
        'independent', {'x'}, ... % y is not actually used in the function but include it so that we can use all 3D data to fit (not just the x and z coordinates)
        'coefficients', {'x1','x2', 'x3'}); % x1 is step location, x2 is pre-step, x3 is post-step

    for aa = 1:length(angles) % rotate x-axis iteratively and fit each time to find best fit
        theta = deg2rad(angles(aa));
        xs_new = xs * cos(theta) - ys * sin(theta); % rotate xs
        ys_new = xs * sin(theta) + ys * cos(theta); % rotate ys

        % Get slice bins
        bins = min(ys_new):slice_size:max(ys_new);
        n_slices = length(bins) - 1;
        rsqrs = nan(n_slices,1);
        for ss = 1:n_slices
            % Get slice data
            slice_bin = bins(ss:ss+1);
            slice_mask = ys_new>=slice_bin(1) & ys_new<slice_bin(2);
            slice_xs = xs_new(slice_mask);
            slice_ys = ys_new(slice_mask);
            slice_Ts = Ts(slice_mask);
            slice_group_data = group_data(slice_mask);
            if length(slice_Ts) < 20
                %disp(['Fewer than 20 data points in slice ' num2str(ss) '/' num2str(n_slices) '... skipping slice']);
                continue
            end

            startpoint_guesses = [mean(slice_xs), mean([slice_group_data(slice_group_data<0);0]), mean([slice_group_data(slice_group_data>0); 0])]; % estimates for step location, pre-step z value, and post-step z value
            [~, gof_curr, ~] = fit(slice_xs, slice_Ts, model, 'StartPoint', startpoint_guesses);
            rsqrs(ss) = gof_curr.rsquare;
        end
        if best_rsquare<mean(rsqrs, 'omitnan', "Weights", metrics_step.numobs)
            best_rsquare = mean(rsqrs, 'omitnan', "Weights", metrics_step.numobs);
            winning_angle = angles(aa);
        end
    end

elseif strcmp(type, 'linear')
    error('linear model not yet implemented in this function')
elseif strcmp(type, 'hinge')
    error('hinge model not yet implemented in this function')
else
    error('model type not recognized')
end
function [winning_angle] = individual_find_axis_noslice(type, xs, ys, Ts, group_data, angles)
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

        startpoint_guesses = [mean(xs_new), mean([group_data(group_data<0);0]), mean([group_data(group_data>0); 0])]; % estimates for step location, pre-step z value, and post-step z value
        [~, gof_curr, ~] = fit(xs_new, Ts, model, 'StartPoint', startpoint_guesses);

        if best_rsquare<gof_curr.rsquare
            best_rsquare = gof_curr.rsquare;
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
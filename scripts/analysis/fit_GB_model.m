function [gof, info, winning_angle, xs_out, ys_out, model_out] = fit_GB_model(type, xs, ys, Ts, group_data, angles)
% FIT_GB_MODEL fits boundary (step function) or gradient (linear, hinge
% functions) to the given data at the supplied angle rotations
%
%   Note that for all custom functions used with Matlab's "fit" function,
%   the default optimizer algorithm is non-linear trust-region

    arguments (Input)
        type {mustBeText, mustBeNonempty} % which model to fit (step, linear, hinge)
        xs (:,1) {mustBeMatrix, mustBeNonNan, mustBeNonempty} % x coordinates
        ys (:,1) {mustBeMatrix, mustBeNonNan, mustBeNonempty} % y coordinates
        Ts (:,1) {mustBeMatrix, mustBeNonNan, mustBeNonempty} % individual data values corresponding to x and y coordinates
        group_data (:,1) {mustBeMatrix, mustBeNonNan, mustBeNonempty} % group data values corresponding to x and y coordinates (used for guessing starting fit values)
        angles {mustBeMatrix, mustBeNonNan, mustBeNonempty} % degrees to rotate the coordinates and fit (this function will output the info for the best fit out of all angles)
    end

    best_rsquare = 0; % initialize
    
    if strcmp(type, 'step')
        model = fittype('x2*(x<x1) + x3*(x>=x1) + y*0',... % step function model
                        'dependent', 'z',...
                        'independent', {'x','y'}, ... % y is not actually used in the function but include it so that we can use all 3D data to fit (not just the x and z coordinates)
                        'coefficients', {'x1','x2', 'x3'}); % x1 is step location, x2 is pre-step, x3 is post-step
        
        for aa = 1:length(angles) % rotate x-axis iteratively and fit each time to find best fit
            theta = deg2rad(angles(aa));
            xs_new = xs * cos(theta) - ys * sin(theta); % rotate xs
            ys_new = xs * sin(theta) + ys * cos(theta); % rotate ys

            startpoint_guesses = [mean(xs_new), mean(group_data(group_data<0)), mean(group_data(group_data>0))]; % estimates for step location, pre-step z value, and post-step z value
            [pre_fit_res, pre_gof, pre_info] = fit([xs_new, ys_new], Ts, model, 'StartPoint', startpoint_guesses);

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
            X = [xs_new, ones(length(xs_new),1)];
            coefs = X \ group_data; % linear regression on only x coordinates
            slope_guess = coefs(1); % 1st param is slope
            intercept_guess = coefs(2); % 2nd is intercept
            startpoint_guesses = [slope_guess, intercept_guess];
    
            [pre_fit_res, pre_gof, pre_info] = fit([xs_new, ys_new], Ts, model, ...
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
            startpoint_guesses = [max([mean(group_data(group_data<0)), min(group_data)]), min([mean(group_data(group_data>0)), max(group_data)]), prctile(xs_new, 25), prctile(xs_new, 75)]; % estimates for
            lower_bounds = [-100, -100, min(xs_new), min(xs_new)];
            upper_bounds = [100, 100, max(xs_new), max(xs_new)];
    
            [pre_fit_res, pre_gof, pre_info] = fit([xs_new, ys_new], Ts, model, ...
                'StartPoint', startpoint_guesses,...
                'Lower', lower_bounds,...
                'Upper', upper_bounds); % fit model
    
            if pre_fit_res.x2< pre_fit_res.x1 % Invalid fit, cannot start the hinge before you end it, retry with bounds that force x1<x2
                [pre_fit_res, pre_gof, pre_info] = fit([xs_new, ys_new], Ts, model, ...
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
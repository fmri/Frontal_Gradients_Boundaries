%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to create boundary and gradient 1D data
%%% with variable gradient and ROI lengths, and fit hinge models to the
%%% data to discern the effect of ROI length on modeled gradient length.
%%%
%%% Tom Possidente - September 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables

% Simulation parameters to explore over
gradient_lengths = 1:1:10; %mm
ROI_lengths = 3:1:20; % total x-axis distance of cortex (mm)

mean_noise_level = 0.5; % will be multiplied with standard normal distribution to determine noise
mean_y_diffs = 3.16; % difference in T-statistics from min to max (derived from average of empirical data)

n_points = 100; % reasonable number of points per slice
N_sims = 100; % number of simulations for each parameter combination

plots = 1;
errors = nan(length(gradient_lengths), length(ROI_lengths));
result = nan(length(gradient_lengths), length(ROI_lengths));

%% Loop over parameter space, simulate data, fit step/hinge models, compare BIC between models, record gradient length error
for rr = 1:length(ROI_lengths)
    ROI_length = ROI_lengths(rr);
    x = linspace(0,ROI_length,n_points);

    parfor gg = 1:length(gradient_lengths)
        gradient_length = gradient_lengths(gg);
        x1 = (ROI_length/2) - (gradient_length/2); % hinge start points
        x2 = (ROI_length/2) + (gradient_length/2); % hinge end point

        metrics_step = nan(N_sims,3); % storage for step fit metrics
        metrics_hinge = nan(N_sims,7); % storage for hinge fit metrics

        for ss = 1:N_sims

            y_diff = mean_y_diffs + randn(); % add some variability to y_diff
            a = -y_diff/2; % min T-stat
            b = y_diff/2; % max T-stat

            noise_level = mean_noise_level + randn()*0.5;
            if noise_level<0.1 % make sure there is some noise at each point
                noise_level = 0.1;
            end
            x_jitter = x' + randn(n_points,1)*0.2; % jitter location of each x coordinate so there are different x coords each time
            y = a.*(x_jitter < x1) + b.*(x_jitter > x2) + ( a + ( (b-a)/(x2-x1)  .* (x_jitter-x1) ) ) .* ( (x_jitter > x1) & (x_jitter < x2) ); % Create underlying data with specified hinge length
            y_noise = y + randn(n_points,1)*noise_level; % add standard normal noise scaled by noise_level
            if plots == 2
                figure;
                scatter(x_jitter,y);
                hold on;
                scatter(x_jitter,y_noise,'filled');
            end

            %% Fit Step model
            model = fittype('x2*(x<x1) + x3*(x>=x1)',... % step function model
                'dependent', 'z',...
                'independent', {'x'}, ...
                'coefficients', {'x1','x2', 'x3'}); % x1 is step location, x2 is pre-step, x3 is post-step

            startpoint_guesses = [mean(x_jitter), -1.5, 1.5]; % if no group data, hardcoded guesses
            [fit_res_step, gof_step, ~] = fit(x_jitter, y_noise, model, 'StartPoint', startpoint_guesses); % fit data
            metrics_step(ss,:) = [gof_step.rsquare, gof_step.rmse, gof_step.sse]; % extract metrics

            %% Fit Hinge model
            model = fittype('a*(x < x1) + b*(x > (x1+x2)) + ( a + ( (b-a)/((x1+x2)-x1)  * (x-x1) ) ) * ( (x > x1) & (x < (x1+x2)) )', ... % piecewise function, constant from -Inf to x1, linear from x1 to x2, constant from x2 to Inf.
                'dependent', 'z',...
                'independent', {'x'}, ...
                'coefficients', {'a','b','x1', 'x2'}); % a is z from -Inf to x1, b is z from x2 to Inf, x1 is where to start linear piece, x2 is where to end linear piece

            startpoint_guesses = [-1.5, 1.5, prctile(x_jitter, 25), prctile(x_jitter, 50)]; % if no group data, hardcoded guesses
            lower_bounds = [-10, -10, min(x_jitter), 0];
            upper_bounds = [10, 10, max(x_jitter), max(x_jitter)-min(x_jitter)];

            [fit_res_hinge, gof_hinge, ~] = fit(x_jitter, y_noise, model, 'StartPoint', startpoint_guesses,...
                'Lower', lower_bounds,...
                'Upper', upper_bounds);  % Fit data
            metrics_hinge(ss,:) = [gof_hinge.rsquare, gof_hinge.rmse, gof_hinge.sse, fit_res_hinge.x1, fit_res_hinge.x2, fit_res_hinge.a, fit_res_hinge.b];

            if plots == 2
                plot([min(x_jitter), metrics_hinge(ss,4), metrics_hinge(ss,4)+metrics_hinge(ss,5), max(x_jitter)],...
                    [metrics_hinge(ss,6), metrics_hinge(ss,6), metrics_hinge(ss,7), metrics_hinge(ss,7)]);
                hold on;
                plot([min(x_jitter), fit_res_step.x1, fit_res_step.x1, max(x_jitter)], ...
                    [fit_res_step.x2, fit_res_step.x2, fit_res_step.x3, fit_res_step.x3])
                title(['True gradient length=' num2str(gradient_length)])
                keyboard;
                close all
            end
        end

        % take mean of metrics across all simulations ('slices')
        hinge_means = mean(metrics_hinge, 'omitnan'); % don't need to weight by # of points in each slice bc they are all the same
        step_means = mean(metrics_step, 'omitnan');

        % Get log likelihoods for each model
        LL_step = loglikelihood(n_points, step_means(2), step_means(3));
        LL_hinge = loglikelihood(n_points, hinge_means(2), hinge_means(3));

        % Compare BICs between models
        [~, bic] = aicbic([LL_step; LL_hinge], [3; 4]); % compare mean BICs

        errors(gg,rr) = gradient_length - hinge_means(5); % expected minus modeled hinge length

        if bic(1)-bic(2)>=10 && hinge_means(1)>0.1 && hinge_means(5)<=4.8 % Hinge model wins BIC and fit, but is less than 4.8mm, therefore boundary wins
            result(gg,rr) = 0;
        elseif bic(1)-bic(2)<=-10 && step_means(1)>0.1 % To win, boundary model must have BIC difference > 10 and r>.1
            result(gg,rr) = 0;
        elseif bic(1)-bic(2)>=10 && hinge_means(1)>0.1 && hinge_means(5)>4.8 % To win, gradient model must have BIC difference > 10 and r>.1, and hinge length > 4.8
            result(gg,rr) = 1;
        else
            result(gg,rr) = nan; % inconclusive
            errors(gg,rr) = nan;
        end

    end

    disp([num2str(rr) '/' num2str(length(ROI_lengths))]);
end

%% Plot
figure;
h = heatmap(result);
h.XDisplayLabels = string(ROI_lengths);
h.YDisplayLabels = string(gradient_lengths);
xlabel('ROI Length');
ylabel('True Gradient Length');
title(['Model Result 1=Gradient, 0=Boundary | Mean Noise=' num2str(mean_noise_level) ' T-stat per point']);
fontsize(16,"points")

figure;
h = heatmap(errors);
h.XDisplayLabels = string(ROI_lengths);
h.YDisplayLabels = string(gradient_lengths);
xlabel('ROI Length');
ylabel('True Gradient Length');
colormap(flipud(redblue_equalized));
clim([-max(abs(errors),[],'all'),max(abs(errors),[],'all')])
title(['Gradient Length Error (mm) | Mean Noise=' num2str(mean_noise_level) ' T-stat per point']);
fontsize(16,"points")

%% Helper functions %%
function LL= loglikelihood(num_obs, rmse, sse)
    LL = -0.5 * num_obs * log( 2*pi*(rmse^2) ) - (1/(2*(rmse^2))) * sse;
end
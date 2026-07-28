%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to create boundary and gradient 1D data,
%%% add varying levels of noise, and fit hinge models to the data to try to
%%% recover the original gradient length.
%%%
%%% Tom Possidente - July 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables

gradient_lengths = 0.01:1:8; %mm
noise_levels = 0:0.4:2;
y_diffs = 2:1:6;

n_points = 100;
x_length = 15;
x = linspace(0,x_length,n_points);
N_sims = 100;

plots = 1;
errors = nan(length(gradient_lengths), length(noise_levels), length(y_diffs));
%%
for gg = 1:length(gradient_lengths)
    gradient_length = gradient_lengths(gg);
    x1 = (x_length/2) - (gradient_length/2);
    x2 = (x_length/2) + (gradient_length/2);
    for yy = 1:length(y_diffs)
        y_diff = y_diffs(yy);
        a = -y_diff/2;
        b = y_diff/2;
        for nn = 1:length(noise_levels)
            noise_level = noise_levels(nn);
            metrics_step = nan(N_sims,3);
            metrics_hinge = nan(N_sims,7);
            for ss = 1:N_sims
                x_jitter = x' + randn(n_points,1)*0.2; % jitter location of each x coordinate so there are different x coords each time
                y = a.*(x_jitter < x1) + b.*(x_jitter > x2) + ( a + ( (b-a)/(x2-x1)  .* (x_jitter-x1) ) ) .* ( (x_jitter > x1) & (x_jitter < x2) );
                y_noise = y + randn(n_points,1)*noise_level;
                if plots == 2
                    figure;
                    scatter(x_jitter,y);
                    hold on;
                    scatter(x_jitter,y_noise,'filled');
                end
                %% Step model
                model = fittype('x2*(x<x1) + x3*(x>=x1)',... % step function model
                    'dependent', 'z',...
                    'independent', {'x'}, ... 
                    'coefficients', {'x1','x2', 'x3'}); % x1 is step location, x2 is pre-step, x3 is post-step

                startpoint_guesses = [mean(x_jitter), -1.5, 1.5]; % if no group data, hardcoded guesses
                [fit_res_step, gof_step, ~] = fit(x_jitter, y_noise, model, 'StartPoint', startpoint_guesses); % actually fit
                metrics_step(ss,:) = [gof_step.rsquare, gof_step.rmse, gof_step.sse];

                %% Hinge model
                model = fittype('a*(x < x1) + b*(x > (x1+x2)) + ( a + ( (b-a)/((x1+x2)-x1)  * (x-x1) ) ) * ( (x > x1) & (x < (x1+x2)) )', ... % piecewise function, constant from -Inf to x1, linear from x1 to x2, constant from x2 to Inf. Must include y in the function even if it has no mathematical effect
                    'dependent', 'z',...
                    'independent', {'x'}, ... 
                    'coefficients', {'a','b','x1', 'x2'}); % a is z from -Inf to x1, b is z from x2 to Inf, x1 is where to start linear piece, x2 is where to end linear piece

                % Calculate initial guesses for gradient model parameters using group-level data
                startpoint_guesses = [-1.5, 1.5, prctile(x_jitter, 25), prctile(x_jitter, 50)]; % if no group data, hardcoded guesses

                lower_bounds = [-10, -10, min(x_jitter), 0];
                upper_bounds = [10, 10, max(x_jitter), max(x_jitter)-min(x_jitter)];

                [fit_res_hinge, gof_hinge, ~] = fit(x_jitter, y_noise, model, 'StartPoint', startpoint_guesses,...
                    'Lower', lower_bounds,...
                    'Upper', upper_bounds); % actually fit
                metrics_hinge(ss,:) = [gof_hinge.rsquare, gof_hinge.rmse, gof_hinge.sse, fit_res_hinge.x1, fit_res_hinge.x2, fit_res_hinge.a, fit_res_hinge.b];
                if plots == 2
                    plot([min(x_jitter), metrics_hinge(ss,4), metrics_hinge(ss,4)+metrics_hinge(ss,5), max(x_jitter)],...
                        [metrics_hinge(ss,6), metrics_hinge(ss,6), metrics_hinge(ss,7), metrics_hinge(ss,7)]);
                    hold on;
                    plot([min(x_jitter), fit_res_step.x1, fit_res_step.x1, max(x_jitter)], ...
                         [fit_res_step.x2, fit_res_step.x2, fit_res_step.x3, fit_res_step.x3])
                    keyboard;
                    close all
                end
            end
            hinge_means = mean(metrics_hinge, 'omitnan');
            step_means = mean(metrics_step, 'omitnan');
            LL_step = loglikelihood(n_points, step_means(2), step_means(3));
            LL_hinge = loglikelihood(n_points, hinge_means(2), hinge_means(3));
            [~, bic] = aicbic([LL_step; LL_hinge], [3; 4]);
            if bic(1)-bic(2)<10 | hinge_means(1)<0.1
                continue
            end
            errors(gg,nn,yy) = gradient_length - hinge_means(5);
        end
    end
    if plots >= 1
        plotdata = squeeze(errors(gg,:,:));
        figure;
        h=heatmap(plotdata);
        h.XDisplayLabels = string(y_diffs);
        h.YDisplayLabels = string(noise_levels);
        xlabel('T stat difference across hinge/boundary');
        ylabel('Noise (scaling factor for standard normal)');
        title(['Underlying Data Gradient length: ' num2str(gradient_length) 'mm']);
    end
end


%% Helper functions %%
function LL= loglikelihood(num_obs, rmse, sse)
LL = -0.5 * num_obs * log( 2*pi*(rmse^2) ) - (1/(2*(rmse^2))) * sse;
end
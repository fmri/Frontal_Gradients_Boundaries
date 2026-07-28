%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to plot some of the visualizations needed
%%% for a methods figure. 
%%%
%%% Tom Possidente - July 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables
% Load in slice model fits for subj AD lh auditory preSMA
load('/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/AD_lh_aud_slice_models.mat', 'metrics_hinge','metrics_step');
N_slices = height(metrics_hinge);
norm_counts = metrics_hinge.numobs/max(metrics_hinge.numobs); % normalize number of data points per slice from 0-1 for easier plotting
hinge_means = mean(metrics_hinge, 'omitnan', 'Weights', metrics_step.numobs); % take mean of all slice models (weighted by number of data points per slice)

%% Plot all slices hinge models
figure;
for ss = 1:N_slices
    if isnan(metrics_hinge{ss,1}) % some slices will be empty/not enough data points
        continue
    else
        p(ss) = plot([metrics_hinge.xmin(ss), metrics_hinge.p1(ss), metrics_hinge.p1(ss)+metrics_hinge.p2(ss), metrics_hinge.xmax(ss)],...
             [metrics_hinge.p3(ss), metrics_hinge.p3(ss), metrics_hinge.p4(ss), metrics_hinge.p4(ss)],...
             'LineWidth', norm_counts(ss)*4, 'Color', [0.8 0.7 0]); % line width corresponds to number of data points in slice
        hold on; 
    end
end

p(ss+1) = plot([hinge_means.xmin, hinge_means.p1, hinge_means.p1+hinge_means.p2, hinge_means.xmax],...
             [hinge_means.p3, hinge_means.p3, hinge_means.p4, hinge_means.p4],...
             'LineWidth', 10, 'Color', 'r', 'LineStyle',':'); % plot mean slice model
ylim([-5,5]);
xlim([65,115]);
xlabel('x');
ylabel('T-stat Difference');
legend([p(9),p(end)], {'Slice Models', 'Subject Mean Model'});
fontsize(16,"points");

%% Plot all subjs mean slice model results

% Load each subject's mean slice model for auditory preSMA LH
load('/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/allsubjs_aud_preSMA_meanslice.mat', 'boundary_wins','good_winning_direction','gradient_wins', 'hinge_means_all');
N_subjs = size(hinge_means_all,1);
hinge_means = [];
figure;
for ss = 1:N_subjs

    if isnan(good_winning_direction(ss,1))
        continue;
    end
    
    % color by gradient/boundary win or inconclusive. 
    if ss==3 % this is subj AD, we want to make it a dashed line to show it is the same subj from the previous plot
        color = [1 0 0];
        style = ':';
    elseif boundary_wins(ss,1)
        color = [0 0.2 1];
        style = '-';
    elseif gradient_wins(ss,1)
        color = [1 0 0];
        style = '-';
    else
        color = [0.7 0.7 0.7];
        style = '-';
    end

    hinge_mean = hinge_means_all{ss,1};
    hinge_means(end+1,:) = table2array(hinge_mean);
    x_modifier = hinge_mean.p1 + (hinge_mean.p2/2); % this will be used to align all subjs on the x-axis to their middle points
    p(ss) = plot([hinge_mean.xmin-x_modifier, hinge_mean.p1-x_modifier, hinge_mean.p1+hinge_mean.p2-x_modifier, hinge_mean.xmax-x_modifier],...
         [hinge_mean.p3, hinge_mean.p3, hinge_mean.p4, hinge_mean.p4], 'Color',color, 'LineWidth',3, 'LineStyle',style);
    hold on; 
end

% Plot group mean model as well
hinge_means = mean(hinge_means);
x_modifier_mean = hinge_means(5) + (hinge_means(6)/2);
p(end+1) = plot([hinge_means(:,9)-x_modifier_mean, hinge_means(:,5)-x_modifier_mean, hinge_means(:,5)+hinge_means(:,6)-x_modifier_mean, hinge_means(10)-x_modifier_mean],...
         [hinge_means(7), hinge_means(7), hinge_means(8), hinge_means(8)], 'Color','k', 'LineWidth',7, "LineStyle", ':');

legend([p(2) p(10) p(5), p(end)], {'Gradient Win', 'Boundary Win', 'Inconclusive', 'Group Mean Model'});
xlabel('x');
ylabel('T-stat Difference');
xlim([-20 20]);
fontsize(16,"points");
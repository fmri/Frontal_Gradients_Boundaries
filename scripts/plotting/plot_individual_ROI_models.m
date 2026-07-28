%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to visualize the single subject winning
%%% gradient/boundary models along the axis of greatest difference of an ROI for each subj
%%%
%%% Tom Possidente - July 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables
ROI = 'preSMA'; %preSMA, inf_lat_frontal, aINS, midIFS, sup_lat_frontal,
modality = 'visual'; % auditory, visual, supramodal
load(['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/Individual_ROI_models/' ROI '_' modality '_hinge_boundary_models.mat'],...
    'boundary_wins', 'gradient_wins', 'hinge_means_all', 'step_means_all', 'good_winning_direction');
N_subjs = size(good_winning_direction,1);
hemis = {'lh','rh'};
N_hemis = length(hemis);

save_png = false;

%% Loop over subjs/hemis and calculate movmean
for hh = 1:N_hemis
    hinge_means = [];
    figure; 
    hemi = hemis{hh};
    for ss = 1:N_subjs

        if isnan(good_winning_direction(ss,hh))
            continue;
        end

        if boundary_wins(ss,hh)
            color = [0 0.2 1];
            style = '-';
        elseif gradient_wins(ss,hh)
            color = [1 0 0];
            style = '-';
        else
            color = [0.7 0.7 0.7];
            style = '-';
        end

        hinge_mean = hinge_means_all{ss,hh};
        hinge_means(end+1,:) = table2array(hinge_mean);
        x_modifier = hinge_mean.p1 + (hinge_mean.p2/2);
        p(ss) = plot([hinge_mean.xmin-x_modifier, hinge_mean.p1-x_modifier, hinge_mean.p1+hinge_mean.p2-x_modifier, hinge_mean.xmax-x_modifier],...
            [hinge_mean.p3, hinge_mean.p3, hinge_mean.p4, hinge_mean.p4], 'Color',color, 'LineWidth',3, 'LineStyle',style);
        hold on;
    end

    hinge_means = mean(hinge_means);
    x_modifier_mean = hinge_means(5) + (hinge_means(6)/2);
    p(end+1) = plot([hinge_means(:,9)-x_modifier_mean, hinge_means(:,5)-x_modifier_mean, hinge_means(:,5)+hinge_means(:,6)-x_modifier_mean, hinge_means(10)-x_modifier_mean],...
        [hinge_means(7), hinge_means(7), hinge_means(8), hinge_means(8)], 'Color','k', 'LineWidth',7, "LineStyle", ':');

    xlabel('x');
    ylabel('T-stat Difference');
    xlim([-20 20]);
    fontsize(16,"points");
    title([ROI ' ' hemi ' ' modality])
    
    if save_png
        print('-dpng', '-r600', [ROI '_' hemi '_' modality '_indsubj_models.png']);
    end
end



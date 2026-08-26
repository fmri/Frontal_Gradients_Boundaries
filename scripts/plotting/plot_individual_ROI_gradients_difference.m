%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to visualize the T-stats diffference 
%%% between WM and sensory contrasts along the axis of greatest difference 
%%% of an ROI for each subj
%%%
%%% Tom Possidente - June 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables
ROI = 'preSMA'; %preSMA, inf_lat_frontal, aINS, midIFS, sup_lat_frontal, aIPS
modality = 'visual'; % auditory, visual, supramodal
load(['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/individual_ROI_tstats/' ROI '_' modality '_WMSMC_xs_Ts.mat'],...
    'Ts_store', 'subjCodes', 'hemis', 'xs_store', 'boundary_x', 'good_winning_direction', 'gradient_wins', 'boundary_wins');
N_subjs = length(subjCodes);
N_hemis = length(hemis);
contrasts = {'sensory', 'WM'};
N_contrasts = 2;

colors = nan(N_subjs,N_hemis,3);
move_mean_window = 2.2; % bin size for averaging (mm)

save_out_image = false;

%% Loop over subjs/hemis and calculate moving mean across ROI
bin_centers = nan(N_subjs, N_hemis, 50);
move_means = nan(N_subjs, N_hemis, 50);
slopes = nan(N_subjs, N_hemis); % also interested in looking at slopes of each

for hh = 1:N_hemis
    hemi = hemis{hh};
    for ss = 1:N_subjs
        subj = subjCodes{ss};
        boundary = boundary_x(ss,hh); % x-axis coordinate (used to align all subjs)
        if gradient_wins(ss,hh) % color by gradient/boundary/inconclusive
            colors(ss,hh,:) = [1 0 0];
        elseif boundary_wins(ss,hh)
            colors(ss,hh,:) = [0 0.2 1];
        elseif good_winning_direction(ss,hh)==0
            colors(ss,hh,:) = [0.7 0.7 0.7];
        else
            colors(ss,hh,:) = [0.7 0.7 0.7];            
        end
        xs = xs_store{ss,hh};
        xs = xs-boundary; % align to boundary point of boundary model on x-axis
        if isempty(xs) % some subjs are missing some ROIs
            continue;
        end
        bin_edges = min(xs):move_mean_window:max(xs);
        bin_centers(ss,hh,1:length(bin_edges)-1) = bin_edges(1:end-1) + (diff(bin_edges)/2);
        Ts_WM = Ts_store{ss,2,hh}; % extract WM T-stats
        Ts_WM(Ts_WM<0) = 0; % clip Tstats at 0 for interpretability of difference T-stats
        Ts_SMC = Ts_store{ss,1,hh};
        Ts_SMC(Ts_SMC<0) = 0;
        Ts_diff = Ts_WM - Ts_SMC;
        coefs = [xs, ones(length(xs),1)] \ Ts_diff'; % linear regression on only x coordinates
        slopes(ss,hh) = coefs(1); % 1st param is slope
        move_mean = nan(length(bin_edges)-1,1);
        % Actually calculate the mean difference T-stat for each bin
        for bb = 1:length(bin_edges)-1
            move_mean(bb) = mean(Ts_diff( (xs>=bin_edges(bb)) & (xs<bin_edges(bb+1)) ));
        end
        move_means(ss,hh, 1:length(move_mean)) = move_mean;
    end
end

%save([ROI '_' modality '_slopes.mat'], 'ROI', 'modality', 'slopes', 'hemis', 'subjCodes', 'contrasts');

%% Get mean of all subjs
means = squeeze(mean(move_means,1,'omitnan'));

for hh = 1:N_hemis
    notnan = sum(~isnan(squeeze(move_means(:,hh,:))), 1);
    means(hh,notnan<5) = nan; % if there are less than 5 subjs with data at a point, make that point nan
end

slope_means = squeeze(mean(slopes,1,'omitnan'));

%% Plot T-stat diffs over x-axis for all subjs (and mean of subjs) for each hemisphere

mean_xs = -24.5*move_mean_window:move_mean_window:24.5*move_mean_window;
mean_xs = squeeze(mean(bin_centers, 1, 'omitnan'))';

f = figure; 
subplot(1,2,1);
for pp = 1:N_subjs
    plot(squeeze(bin_centers(pp,1,:))', squeeze(move_means(pp,1,:))', '-', 'Color', squeeze(colors(pp,1,:)), 'LineWidth',1.25);
    hold on; 
end
plot(mean_xs(:,1), squeeze(means(1,:)), '-*', 'Color', 'k', 'LineWidth',4);
title('LH');
ylabel('WM-Sensory T-stats');
ylim([min(move_means(:,1,:),[],'all')-0.2, max(move_means(:,1,:),[],'all')+0.2])
xlim([min(bin_centers(:,1,:),[],'all'), max(bin_centers(:,1,:),[],'all')])

subplot(1,2,2);
for pp = 1:N_subjs
    plot(squeeze(bin_centers(pp,2,:))', squeeze(move_means(pp,2,:))', '-', 'Color', squeeze(colors(pp,2,:)), 'LineWidth',1.25);
    hold on; 
end
plot(mean_xs(:,2), squeeze(means(2,:)), '-*', 'Color', 'k', 'LineWidth',4);
title('RH');
ylim([min(move_means(:,2,:),[],'all')-0.2, max(move_means(:,2,:),[],'all')+0.2])
xlim([min(bin_centers(:,2,:),[],'all'), max(bin_centers(:,2,:),[],'all')])

fontsize(14,"points")
sgtitle([replace(ROI,'_', ' ') ' | ' modality]);
if save_out_image
    print('-dpng', '-r600', ['tstat_diff_lines_' ROI '_' modality '.png'])
end
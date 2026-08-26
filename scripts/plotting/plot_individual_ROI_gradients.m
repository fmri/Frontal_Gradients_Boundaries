%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to visualize the WM and sensory contrast
%%% T-stats along the axis of greatest difference of an ROI for each subj
%%%
%%% Tom Possidente - May 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables
ROI = 'preSMA'; %preSMA, inf_lat_frontal, aINS, midIFS, sup_lat_frontal, 
modality = 'visual'; % auditory, visual, supramodal
load(['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/individual_ROI_tstats/' ROI '_' modality '_WMSMC_xs_Ts.mat'],...
    'Ts_store', 'subjCodes', 'hemis', 'xs_store', 'boundary_x', 'good_winning_direction', 'gradient_wins', 'boundary_wins');
N_subjs = length(subjCodes);
N_hemis = length(hemis);
contrasts = {'sensory', 'WM'};
N_contrasts = 2;

move_mean_window = 2.2; % bin size for averaging (mm)

save_out_image = true;

%% Loop over subjs/hemis and calculate moving mean across ROI
bin_centers = nan(N_subjs, N_hemis, 50);
move_means = nan(N_subjs, N_contrasts, N_hemis, 50);
slopes = nan(N_subjs, N_contrasts, N_hemis); % also interested in looking at slopes of each

for hh = 1:N_hemis
    hemi = hemis{hh};
    for ss = 1:N_subjs
        subj = subjCodes{ss};
        boundary = boundary_x(ss,hh); % x-axis coordinate (used to align all subjs)
        xs = xs_store{ss,hh};
        xs = xs-boundary; % align to boundary point of boundary model on x-axis
        if isempty(xs) 
            continue;
        end
        bin_edges = min(xs):move_mean_window:max(xs);
        bin_centers(ss,hh,1:length(bin_edges)-1) = bin_edges(1:end-1) + (diff(bin_edges)/2);
        for cc = 1:N_contrasts
            Ts = Ts_store{ss,cc,hh};
            coefs = [xs, ones(length(xs),1)] \ Ts'; % linear regression on only x coordinates
            slopes(ss,cc,hh) = coefs(1); % 1st param is slope
            move_mean = nan(length(bin_edges)-1,1);
            % Actually calculate the mean difference T-stat for each bin
            for bb = 1:length(bin_edges)-1
                move_mean(bb) = mean(Ts( (xs>=bin_edges(bb)) & (xs<bin_edges(bb+1)) ));
            end
            move_means(ss,cc,hh, 1:length(move_mean)) = move_mean;
        end
    end
end

%save([ROI '_' modality '_slopes.mat'], 'ROI', 'modality', 'slopes', 'hemis', 'subjCodes', 'contrasts');

%% Get means 
means = squeeze(mean(move_means,1,'omitnan'));

for hh = 1:N_hemis
    for cc = 1:N_contrasts
        notnan = sum(~isnan(squeeze(move_means(:,cc,hh,:))), 1);
        means(cc,hh,notnan<5) = nan; % if there are less than 5 subjs with data at a point, make that point nan
    end
end

slope_means = squeeze(mean(slopes,1,'omitnan'));

%% Plot T-stats for WM and Sensory (SMC) over x-axis for all subjs (and mean of subjs) for each hemisphere

mean_xs = squeeze(mean(bin_centers, 1, 'omitnan'))';

f = figure; 
subplot(3,2,1);
plot(squeeze(bin_centers(:,1,:))', squeeze(move_means(:,1,1,:))', '-', 'Color', '#32CD32');
hold on; 
plot(mean_xs(:,1), squeeze(means(1,1,:)), '-*', 'Color', 'k', 'LineWidth',4);
title('LH');
ylabel('Sensory Drive T-stats');
ylim([min(move_means(:,1,1,:),[],'all')-0.2, max(move_means(:,1,1,:),[],'all')+0.2])
xlim([min(bin_centers(:,1,:),[],'all'), max(bin_centers(:,1,:),[],'all')])
subplot(3,2,5);
p1 = plot(squeeze(bin_centers(:,1,:))', squeeze(move_means(:,1,1,:))', '-', 'Color', '#32CD32');
hold on; 
plot(mean_xs(:,1), squeeze(means(1,1,:)), '-*', 'Color', 'k', 'LineWidth',4)
ylabel('Both');
ylim([min( [min(move_means(:,1,1,:),[],'all'), min(move_means(:,2,1,:),[],'all')] )-0.2, ...
      max( [max(move_means(:,1,1,:),[],'all'), max(move_means(:,2,1,:),[],'all')] )+0.2 ])
xlim([min(bin_centers(:,1,:),[],'all'), max(bin_centers(:,1,:),[],'all')])

subplot(3,2,2);
plot(squeeze(bin_centers(:,2,:))', squeeze(move_means(:,1,2,:))', '-', 'Color', '#32CD32');
hold on; 
plot(mean_xs(:,2), squeeze(means(1,2,:)), '-*', 'Color', 'k', 'LineWidth',4);
title('RH');
ylim([min(move_means(:,1,2,:),[],'all')-0.2, max(move_means(:,1,2,:),[],'all')+0.2])
xlim([min(bin_centers(:,2,:),[],'all'), max(bin_centers(:,2,:),[],'all')])
subplot(3,2,6);
plot(squeeze(bin_centers(:,2,:))', squeeze(move_means(:,1,2,:))', '-', 'Color', '#32CD32');
hold on; 
ylim([min( [min(move_means(:,1,2,:),[],'all'), min(move_means(:,2,2,:),[],'all')] )-0.2, ...
      max( [max(move_means(:,1,2,:),[],'all'), max(move_means(:,2,2,:),[],'all')] )+0.2 ])
xlim([min(bin_centers(:,2,:),[],'all'), max(bin_centers(:,2,:),[],'all')])

subplot(3,2,3);
plot(squeeze(bin_centers(:,1,:))', squeeze(move_means(:,2,1,:))', '-', 'Color', 'r');
hold on; 
plot(mean_xs(:,1), squeeze(means(2,1,:)), '-*', 'Color', [0.4,0.4,0.4], 'LineWidth',4);
ylabel('WM T-stats')
ylim([min(move_means(:,2,1,:),[],'all')-0.2, max(move_means(:,2,1,:),[],'all')+0.2])
xlim([min(bin_centers(:,1,:),[],'all'), max(bin_centers(:,1,:),[],'all')])
subplot(3,2,5);
p2 = plot(squeeze(bin_centers(:,1,:))', squeeze(move_means(:,2,1,:))', '-', 'Color', 'r');
hold on; 
p3=plot(mean_xs(:,1), squeeze(means(1,1,:)), '-*', 'Color', [0,0,0], 'LineWidth',4);
p4=plot(mean_xs(:,1), squeeze(means(2,1,:)), '-*', 'Color', [0.4,0.4,0.4], 'LineWidth',4);
xlabel('Cortical Space (mm)');
ylim([min( [min(move_means(:,2,1,:),[],'all'), min(move_means(:,1,1,:),[],'all')] )-0.2, ...
      max( [max(move_means(:,2,1,:),[],'all'), max(move_means(:,1,1,:),[],'all')] )+0.2 ])
xlim([min(bin_centers(:,1,:),[],'all'), max(bin_centers(:,1,:),[],'all')])
legend([p1(1),p2(1),p3,p4], {'individual sensory', 'individual WM' 'mean sensory', 'mean WM'}, 'Position',[0.42682007778216,0.004234291861436,0.169856459330144,0.072005383580081]);

subplot(3,2,4);
plot(squeeze(bin_centers(:,2,:))', squeeze(move_means(:,2,2,:))', '-', 'Color', 'r');
hold on; 
plot(mean_xs(:,2), squeeze(means(2,2,:)), '-*', 'Color', [0.4,0.4,0.4], 'LineWidth',4)
ylim([min(move_means(:,2,2,:),[],'all')-0.2, max(move_means(:,2,2,:),[],'all')+0.2])
xlim([min(bin_centers(:,2,:),[],'all'), max(bin_centers(:,2,:),[],'all')])
subplot(3,2,6);
plot(squeeze(bin_centers(:,2,:))', squeeze(move_means(:,2,2,:))', '-', 'Color', 'r');
hold on; 
plot(mean_xs(:,2), squeeze(means(1,2,:)), '-*', 'Color', 'k', 'LineWidth',4);
plot(mean_xs(:,2), squeeze(means(2,2,:)), '-*', 'Color', [0.4,0.4,0.4], 'LineWidth',4);
xlabel('Cortical Space (mm)');
ylim([min( [min(move_means(:,1,2,:),[],'all'), min(move_means(:,2,2,:),[],'all')] )-0.2, ...
      max( [max(move_means(:,1,2,:),[],'all'), max(move_means(:,2,2,:),[],'all')] )+0.2 ])

sgtitle([replace(ROI,'_', ' ') ' | ' modality]);
set(f,"Position",[1 1 782 1000]);
if save_out_image
    print('-dpng', '-r600', ['tstat_lines_' ROI '_' modality '.png'])
end
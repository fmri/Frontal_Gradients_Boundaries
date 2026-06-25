%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to visualize the T-stats along the axis of
%%% greatest difference of an ROI for each subj
%%%
%%% Tom Possidente - June 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables
ROI = 'preSMA'; %preSMA, inf_lat_frontal, aINS, midIFS, sup_lat_frontal, aIPS
modality = 'auditory'; % auditory, visual, supramodal
load(['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/individual_ROI_tstats/' ROI '_' modality '_WMSMC_xs_Ts.mat'],...
    'Ts_store', 'subjCodes', 'hemis', 'xs_store', 'boundary_x');
N_subjs = length(subjCodes);
N_hemis = length(hemis);
contrasts = {'sensory', 'WM'};
N_contrasts = 2;

move_mean_window = 2.2; % mm

save_out_image = true;

%% Loop over subjs/hemis and calculate movmean
bin_centers = nan(N_subjs, N_hemis, 50);
move_means = nan(N_subjs, N_hemis, 50);
slopes = nan(N_subjs, N_hemis);

for hh = 1:N_hemis
    hemi = hemis{hh};
    for ss = 1:N_subjs
        subj = subjCodes{ss};
        boundary = boundary_x(ss,hh);
        xs = xs_store{ss,hh};
        xs = xs-boundary; % align to start at 0
        if isempty(xs)
            continue;
        end
        bin_edges = min(xs):move_mean_window:max(xs);
        bin_centers(ss,hh,1:length(bin_edges)-1) = bin_edges(1:end-1) + (diff(bin_edges)/2);
        Ts_WM = Ts_store{ss,1,hh};
        Ts_WM(Ts_WM<0) = 0;
        Ts_SMC = Ts_store{ss,2,hh};
        Ts_SMC(Ts_WM<0) = 0;
        Ts_diff = Ts_WM - Ts_SMC;
        coefs = [xs, ones(length(xs),1)] \ Ts_diff'; % linear regression on only x coordinates
        slopes(ss,hh) = coefs(1); % 1st param is slope
        move_mean = nan(length(bin_edges)-1,1);
        for bb = 1:length(bin_edges)-1
            move_mean(bb) = mean(Ts_diff( (xs>=bin_edges(bb)) & (xs<bin_edges(bb+1)) ));
        end
        move_means(ss,hh, 1:length(move_mean)) = move_mean;
    end
end

%save([ROI '_' modality '_slopes.mat'], 'ROI', 'modality', 'slopes', 'hemis', 'subjCodes', 'contrasts');

%% Get means 
means = squeeze(mean(move_means,1,'omitnan'));

for hh = 1:N_hemis
    notnan = sum(~isnan(squeeze(move_means(:,hh,:))), 1);
    means(hh,notnan<5) = nan;
end

slope_means = squeeze(mean(slopes,1,'omitnan'));

%% Plot

mean_xs = -24.5*move_mean_window:move_mean_window:24.5*move_mean_window;
mean_xs = squeeze(mean(bin_centers, 1, 'omitnan'))';

f = figure; 
subplot(1,2,1);
plot(squeeze(bin_centers(:,1,:))', squeeze(move_means(:,1,:))', '-', 'Color', 'b');
hold on; 
plot(mean_xs(:,1), squeeze(means(1,:)), '-*', 'Color', 'k', 'LineWidth',4);
title('LH');
ylabel('WM-Sensory T-stats');
ylim([min(move_means(:,1,:),[],'all')-0.2, max(move_means(:,1,:),[],'all')+0.2])
xlim([min(bin_centers(:,1,:),[],'all'), max(bin_centers(:,1,:),[],'all')])

subplot(1,2,2);
plot(squeeze(bin_centers(:,2,:))', squeeze(move_means(:,2,:))', '-', 'Color', 'b');
hold on; 
plot(mean_xs(:,2), squeeze(means(2,:)), '-*', 'Color', 'k', 'LineWidth',4);
title('RH');
ylim([min(move_means(:,2,:),[],'all')-0.2, max(move_means(:,2,:),[],'all')+0.2])
xlim([min(bin_centers(:,2,:),[],'all'), max(bin_centers(:,2,:),[],'all')])

fontsize(14,"points")
sgtitle([replace(ROI,'_', ' ') ' | ' modality]);
%set(f,"Position",[1 1 782 1000]);
if save_out_image
    print('-dpng', '-r600', ['tstat_diff_lines_' ROI '_' modality '.png'])
end
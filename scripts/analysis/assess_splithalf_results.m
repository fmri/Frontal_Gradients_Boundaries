%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to plot and run statistical tests on the
%%% split-half cross validation gradient vs boundary model fits
%%%
%%% Tom Possidente - Septemeber 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables

ROIs = {'preSMA', 'inf_lat_frontal', 'aINS', 'sup_lat_frontal', 'midIFS', 'aIPS'};
modalities = {'auditory', 'visual'};

data_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/cross_validation_modelcomparison_results/';

results = table();

%% Loop over ROIs and modalities, extract split-half results
for rr = 1:length(ROIs)
    ROI = ROIs{rr};
    for mm = 1:length(modalities)
        modality = modalities{mm};
        res = load([data_dir 'CV_modelcomp_results_' ROI '_' modality '_nobaddir.mat']);
        linear_xdist_hemimean = mean(res.linear_xdist,2,'omitnan');
        [~,hinge_length_pval] = ttest(linear_xdist_hemimean, 4.8); 
        mean_hingedist = mean(linear_xdist_hemimean, 'omitnan');
        SE_hingedist = std(linear_xdist_hemimean,'omitnan')/sqrt(sum(~isnan(linear_xdist_hemimean)));
        results = [results; {ROI, modality, res.p, hinge_length_pval, res.mean_step_CVcorr, res.mean_hinge_CVcorr, ...
                             res.SE_step_CVcorr, res.SE_hinge_CVcorr, mean_hingedist, SE_hingedist}];
    end
end

results.Properties.VariableNames = {'ROI', 'modality', 'pval_corr', 'pval_length', 'mean_step_CVcorr', 'mean_hinge_CVcorr',...
                                    'SE_step_CVcorr', 'SE_hinge_CVcorr', 'mean_hingedist', 'SE_hingedist'};

%% FDR BH hypothesis testing
[corrected_p_corr, ~, h_corr] = fdr_BH(results.pval_corr, 0.05); % testing CV corr values from gradient vs boundary models
results.corrected_pval_corr = corrected_p_corr';
results.h_corr = h_corr;

[corrected_p_length, ~, h_length] = fdr_BH(results.pval_length, 0.05); % testing CV gradient lengths from hinge models
results.corrected_pval_length = corrected_p_length';
results.h_length = h_length;

%% Plot
figure;
subplot(2,2,1);
vis_corrdata = results(strcmp(results.modality, 'visual'),5:8);
grouped_errorbars(vis_corrdata{:,1:2}, vis_corrdata{:,3:4});
ylabel('Correlation');
%xticklabels(strrep(ROIs,'_',' '));
title('Visual Cross-validation');

subplot(2,2,2);
aud_corrdata = results(strcmp(results.modality, 'auditory'),5:8);
[bar_handles, errbar_handles] = grouped_errorbars(aud_corrdata{:,1:2}, aud_corrdata{:,3:4});
%xticklabels(strrep(ROIs,'_',' '));
title('Auditory Cross-validation');
legend({'Boundary Model','Gradient Model'})

subplot(2,2,3);
vis_hingedist = results.mean_hingedist(strcmp(results.modality,'visual'));
bar(vis_hingedist);
hold on;
errorbar(vis_hingedist, results.SE_hingedist(strcmp(results.modality,'visual')), '.k');
%xticklabels(strrep(ROIs,'_',' '));
ylabel('Gradient Distance (mm)');
title('Visual Gradient Distance');
yline(4.8,'--r','LineWidth',3);
ylim([0,11]);

subplot(2,2,4);
aud_hingedist = results.mean_hingedist(strcmp(results.modality,'auditory'));
bar(aud_hingedist);
hold on;
errorbar(aud_hingedist, results.SE_hingedist(strcmp(results.modality,'auditory')), '.k');
%xticklabels(strrep(ROIs,'_',' '));
ylabel('Gradient Distance (mm)');
title('Auditory Gradient Distance');
yline(4.8,'--r','LineWidth',3);
ylim([0,11]);

fontsize(16,'points');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is run statistical tests on model fit
%%% (r-squared) and hinge length parameters from individual model fitting
%%%
%%% Tom Possidente - June 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables
data_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/r2_hingelength_data/';
ROIs = {'preSMA', 'inf_lat_frontal', 'aINS', 'sup_lat_frontal', 'aIPS', 'midIFS'};
N_ROIs = length(ROIs);
ROI_hemis_vis = {'LR', 'LR', 'L', 'L', '', 'L'};
ROI_hemis_aud = {'LR', 'LR', 'LR', 'L', 'R', 'L'};
ROI_hemis = {ROI_hemis_aud; ROI_hemis_vis};
modalities = {'auditory', 'visual'};
N_modalities = length(modalities);

%% Loop over ROIs and modality and run LMEs
rsqr_estimates = nan(N_ROIs,N_modalities);
rsqr_SEs = nan(N_ROIs,N_modalities);
rsqr_pvals = nan(N_ROIs,N_modalities);
hingedist_estimates = nan(N_ROIs,N_modalities);
hingedist_SEs = nan(N_ROIs,N_modalities);
hingedist_pvals = nan(N_ROIs,N_modalities);
names = cell(N_ROIs,N_modalities);
for rr = 1:N_ROIs
    ROI = ROIs{rr};
    for mm = 1:N_modalities
        modality = modalities{mm};
        names{rr,mm} = [ROI '_' modality];
        if isfile([data_dir ROI '_' modality '_r2_hingelength.mat'])
            load([data_dir ROI '_' modality '_r2_hingelength.mat'], 'gradient_wins', 'boundary_wins', 'good_winning_direction', 'winning_rsquare', 'linear_xdist');
        else
            disp([ROI ' ' modality ' file not found']);
            continue;
        end

        % Testing r-squared > .1 and hinge length > 4.8
        ROI_hemi = ROI_hemis{mm}{rr};
        rsqr_table = table();
        hinge_table = table; 
        if contains(ROI_hemi, 'L')
            % r2
            subj_ids_r2_lh = 1:size(winning_rsquare,1);    
            rsqr_lh = winning_rsquare(gradient_wins(:,1),1);
            subj_ids_r2_lh = subj_ids_r2_lh(gradient_wins(:,1))';
            rsqr_table = table(rsqr_lh, subj_ids_r2_lh, ones(length(rsqr_lh),1), 'VariableNames',{'rsqr', 'subj', 'hemi'});
            % hinge
            subj_ids_hinge_lh = 1:size(linear_xdist,1);
            hinge_lh = linear_xdist(~isnan(linear_xdist(:,1))&good_winning_direction(:,1)==1,1);
            subj_ids_hinge_lh = subj_ids_hinge_lh(~isnan(linear_xdist(:,1))&good_winning_direction(:,1)==1)';
            hinge_table = table(hinge_lh, subj_ids_hinge_lh, ones(length(hinge_lh),1), 'VariableNames',{'hinge', 'subj', 'hemi'});
        end
        if contains(ROI_hemi, 'R')
            % r2
            subj_ids_rh = 1:size(winning_rsquare,1);
            rsqr_rh = winning_rsquare(gradient_wins(:,2),2);
            subj_ids_rh = subj_ids_rh(gradient_wins(:,2))';
            rsqr_table_rh = table(rsqr_rh, subj_ids_rh, 2*ones(length(rsqr_rh),1), 'VariableNames',{'rsqr', 'subj', 'hemi'});
            rsqr_table = [rsqr_table; rsqr_table_rh];
            % hinge
            subj_ids_hinge_rh = 1:size(linear_xdist,1);
            hinge_rh = linear_xdist(~isnan(linear_xdist(:,2))&good_winning_direction(:,2)==1,2);
            subj_ids_hinge_rh = subj_ids_hinge_rh(~isnan(linear_xdist(:,2))&good_winning_direction(:,2)==1)';
            hinge_table_rh = table(hinge_rh, subj_ids_hinge_rh, 2*ones(length(hinge_rh),1), 'VariableNames',{'hinge', 'subj', 'hemi'});
            hinge_table = [hinge_table; hinge_table_rh];
        end

        disp([ROI ' ' modality])

        
        rsqr_table.subj = categorical(rsqr_table.subj);
        rsqr_table.hemi = categorical(rsqr_table.hemi);
        rsqr_table.rsqr = rsqr_table.rsqr - 0.1; % shift down 0.1 so testing against 0 is the same as testing against 0.1 (null hypothesis)
        hinge_table.subj = categorical(hinge_table.subj);
        hinge_table.hemi = categorical(hinge_table.hemi);
        hinge_table.hinge = hinge_table.hinge - 4.8; % shift down 4.8 so testing against 0 is same as testing against 4.8

        if length(ROI_hemi)==2
            rsqr_model = fitlme(rsqr_table, 'rsqr ~ 1 + (1|hemi) + (1|subj)');
            hinge_model = fitlme(hinge_table, 'hinge ~ 1 + (1|hemi) + (1|subj)');
        else
            rsqr_model = fitlme(rsqr_table, 'rsqr ~ 1');
            hinge_model = fitlme(hinge_table, 'hinge ~ 1');
        end

        rsqr_estimates(rr,mm) = rsqr_model.Coefficients.Estimate(1);
        rsqr_SEs(rr,mm) = rsqr_model.Coefficients.SE(1);
        rsqr_pvals(rr,mm) = rsqr_model.Coefficients.pValue(1);
        if length(rsqr_model.Coefficients.pValue) ==2
            if rsqr_model.Coefficients.pValue<0.1
                disp([ROI ' ' modality ' r^2 significant hemisphere difference'])
            end
        end

        hingedist_estimates(rr,mm) = hinge_model.Coefficients.Estimate(1);
        hingedist_SEs(rr,mm) = hinge_model.Coefficients.SE(1);
        hingedist_pvals(rr,mm) = hinge_model.Coefficients.pValue(1);
        if length(hinge_model.Coefficients.pValue) ==2
            if hinge_model.Coefficients.pValue<0.1
                disp([ROI ' ' modality ' hinge significant hemisphere difference'])
            end
        end
    end
end

[~,~,h] = fdr_BH([rsqr_pvals/2; hingedist_pvals/2],0.05,false); % divide each pval by 2 because we have a directional hypotheses and want a one-tailed test whereas fitlme automatically does a 2-tailed test
h = reshape(h, 12,2);
h(1:6,:)
names
h(7:end,:)

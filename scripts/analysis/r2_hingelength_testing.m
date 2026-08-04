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
ROI_hemis_vis = {'C', 'C', 'LR', 'LR', '', 'LR'}; % hemispheres being examined for each ROI (for visual)
ROI_hemis_aud = {'C', 'C', 'C', 'C', 'LR', 'LR'}; % hemispheres being examined for each ROI (for auditory)
ROI_hemis = {ROI_hemis_aud; ROI_hemis_vis};
modalities = {'auditory', 'visual'};
N_modalities = length(modalities);
hemis = {'lh', 'rh'};
N_hemis = length(hemis);

%% Loop over ROIs and modality and run LMEs
pval_table = table();
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
        % LH r2
        subj_ids_r2_lh = 1:size(winning_rsquare,1);
        if sum(gradient_wins(:,1))>sum(boundary_wins(:,1))
            winningest_model_wins = gradient_wins(:,1);
        elseif sum(gradient_wins(:,1))<sum(boundary_wins(:,1))
            keyboard;
            winningest_model_wins = boundary_wins(:,1);
        else
            error('equal number of gradient and boundary wins (I dont think this ever happens)');
        end
        rsqr_lh = winning_rsquare(winningest_model_wins&good_winning_direction(:,1)==1,1);
        subj_ids_r2_lh = subj_ids_r2_lh(winningest_model_wins&good_winning_direction(:,1)==1)';
        rsqr_table = table(rsqr_lh, subj_ids_r2_lh, ones(length(rsqr_lh),1), 'VariableNames',{'rsqr', 'subj', 'hemi'});
        % LH hinge
        subj_ids_hinge_lh = 1:size(linear_xdist,1);
        hinge_lh = linear_xdist(~isnan(linear_xdist(:,1))&good_winning_direction(:,1)==1,1);
        subj_ids_hinge_lh = subj_ids_hinge_lh(~isnan(linear_xdist(:,1))&good_winning_direction(:,1)==1)';
        hinge_table = table(hinge_lh, subj_ids_hinge_lh, ones(length(hinge_lh),1), 'VariableNames',{'hinge', 'subj', 'hemi'});
        % RH r2
        subj_ids_rh = 1:size(winning_rsquare,1);
        if sum(gradient_wins(:,2))>sum(boundary_wins(:,2))
            winningest_model_wins = gradient_wins(:,2);
        elseif sum(gradient_wins(:,2))<sum(boundary_wins(:,2))
            keyboard;
            winningest_model_wins = boundary_wins(:,2);
        else
            error('equal number of gradient and boundary wins (I dont think this ever happens)');
        end
        rsqr_rh = winning_rsquare(winningest_model_wins&good_winning_direction(:,2)==1,2);
        subj_ids_rh = subj_ids_rh(winningest_model_wins&good_winning_direction(:,2)==1)';
        rsqr_table_rh = table(rsqr_rh, subj_ids_rh, 2*ones(length(rsqr_rh),1), 'VariableNames',{'rsqr', 'subj', 'hemi'});
        rsqr_table = [rsqr_table; rsqr_table_rh];
        % RH hinge
        subj_ids_hinge_rh = 1:size(linear_xdist,1);
        hinge_rh = linear_xdist(~isnan(linear_xdist(:,2))&good_winning_direction(:,2)==1,2);
        subj_ids_hinge_rh = subj_ids_hinge_rh(~isnan(linear_xdist(:,2))&good_winning_direction(:,2)==1)';
        hinge_table_rh = table(hinge_rh, subj_ids_hinge_rh, 2*ones(length(hinge_rh),1), 'VariableNames',{'hinge', 'subj', 'hemi'});
        hinge_table = [hinge_table; hinge_table_rh];

        disp([ROI ' ' modality])

        rsqr_table.subj = categorical(rsqr_table.subj); % make these variables categorical instead of continuous
        rsqr_table.hemi = categorical(rsqr_table.hemi);
        rsqr_table.rsqr = rsqr_table.rsqr - 0.1; % shift down 0.1 so testing against 0 is the same as testing against 0.1 (null hypothesis)
        hinge_table.subj = categorical(hinge_table.subj);
        hinge_table.hemi = categorical(hinge_table.hemi);
        hinge_table.hinge = hinge_table.hinge - 4.8; % shift down 4.8 so testing against 0 is same as testing against 4.8

        if length(ROI_hemi)==2 % test hemispheres separately to see if there is a hemispheric difference
            for hh = 1:N_hemis
                rsqr_model = fitlme(rsqr_table(rsqr_table.hemi==string(hh),:), 'rsqr ~ 1');
                hinge_model = fitlme(hinge_table(hinge_table.hemi==string(hh),:), 'hinge ~ 1');
                pval_table = [pval_table; cell2table({'rsqr', hemis{hh}, modality, ROI, rsqr_model.Coefficients.Estimate+0.1, ...
                              rsqr_model.Coefficients.SE, rsqr_model.Coefficients.pValue(1)})];
                pval_table = [pval_table; cell2table({'hinge', hemis{hh}, modality, ROI, hinge_model.Coefficients.Estimate+4.8, ...
                              hinge_model.Coefficients.SE, hinge_model.Coefficients.pValue(1)})];
            end
        else % combine hemispheres and use hemi as a random effect (and subj)
            rsqr_model = fitlme(rsqr_table, 'rsqr ~ 1 + (1|hemi) + (1|subj)');
            hinge_model = fitlme(hinge_table, 'hinge ~ 1 + (1|hemi) + (1|subj)');
            pval_table = [pval_table; cell2table({'rsqr', 'C', modality, ROI, rsqr_model.Coefficients.Estimate+0.1, ...
                              rsqr_model.Coefficients.SE, rsqr_model.Coefficients.pValue(1)})];
            pval_table = [pval_table; cell2table({'hinge', 'C', modality, ROI, hinge_model.Coefficients.Estimate+4.8, ...
                              hinge_model.Coefficients.SE, hinge_model.Coefficients.pValue(1)})];
        end
    end
end
pval_table.Properties.VariableNames = {'metric', 'hemisphere', 'modality', 'ROI', 'estimate', 'SE', 'pval'};
[~,~,h] = fdr_BH(pval_table.pval/2,0.05,false); % divide each pval by 2 because we have a directional hypotheses and want a one-tailed test whereas fitlme automatically does a 2-tailed test
pval_table.sig_fdrbh = h

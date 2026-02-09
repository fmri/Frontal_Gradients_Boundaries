%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to compute sensory drive and working
%%% memory activation separately for each ROI across the cortex.
%%%
%%%
%%% Tom Possidente - January 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables
ROIs = {'aMFG', 'midIFS', 'aINS', 'preSMA', 'inf_lat_frontal', 'sup_lat_frontal', 'midINS', 'ant_temporal', 'parietal_opercular', 'post_temporal', 'aIPS', ...
    'ms_post_STSG', 'VOT', 'cIPS', 'MT', 'LOT', 'pIPS', 'VO', 'DO', 'post_col_sulc'};
N_ROIs = length(ROIs);

hemis = {'lh', 'rh'};

subjCodes = {'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS', 'AI'};
N = length(subjCodes);

contrasts = {'vP-f', 'vA-vP', 'aP-f', 'aA-aP', 'aPvP-f', 'vAaA-vPaP'};
contrast_types = {'visual', 'visual', 'auditory', 'auditory', 'supramodal', 'supramodal'};
N_contrasts = length(contrasts);

N_vertices = 163842;
vertex_inds = 1:N_vertices;

ROI_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/subj_specific_01/';
data_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';


%% Loop through subjs/ROIs/hemi/contrasts and extract mean activation
PSCs = nan(N, N_ROIs, 2, N_contrasts);
coordinates = nan(N, N_ROIs, 2, N_contrasts);

for ss = 1:N
    subjCode = subjCodes{ss};
    for hh = 1:2
        hemi = hemis{hh};
        for cc = 1:N_contrasts
            contrast = contrasts{cc};
            contrast_type = contrast_types{cc};

            if strcmp(contrast, 'aP-f') | strcmp(contrast, 'vP-f') % these contrasts are backwards
                split = strsplit(contrast,'-');
                contrast = [split{2} '-' split{1}];
                reverse_contrast = true;
            else
                reverse_contrast = false;
            end
            
            data = MRIread([data_dir subjCode '/localizer/localizer_contrasts_' hemi '/' contrast '/cespct.nii.gz']);

            if reverse_contrast
                data.vol = -data.vol;
            end

            for rr = 1:N_ROIs
                ROI_name = ROIs{rr};
                label_fpath = [ROI_dir hemi '.' ROI_name '_' contrast_type '_' subjCode '.label'];
                if isfile(label_fpath) % if file doesn't exist, leave PSC as nan
                    ROI_label = readtable(label_fpath, 'FileType','text');
                    coordinates(ss,rr,hh,cc) = mean(ROI_label{:,3});
                    ROI_label_mask = ismember(vertex_inds, ROI_label{:,1}+1);
                    ROI_data = data.vol(ROI_label_mask);
                    PSCs(ss,rr,hh,cc) = mean(ROI_data);
                end
            end
        end
    end
end

%% Check distributions of PSCs
for cc = 1:N_contrasts
    data = PSCs(:,:,:,cc);
    figure; 
    histogram(data(:));
    xlabel('PSCs');
    title(contrasts{cc});
end

%% Threshold negatives to near zero
PSCs(PSCs<0) = eps;

%% Plot Average PSC ratios for each ROI
PSC_ratios = nan(N,N_ROIs, 2, 3);
PSC_ratios(:,:,:,1) = PSCs(:,:,:,2) ./ (PSCs(:,:,:,1)+PSCs(:,:,:,2)); % WM / SMC+WM visual
PSC_ratios(:,:,:,2) = PSCs(:,:,:,4) ./ (PSCs(:,:,:,3)+PSCs(:,:,:,4)); % WM / SMC+WM auditory
PSC_ratios(:,:,:,3) = PSCs(:,:,:,6) ./ (PSCs(:,:,:,5)+PSCs(:,:,:,6)); % WM / SMC+WM supramodal

means = squeeze(mean(PSC_ratios, [1,3], 'omitnan'));
means_tbl = table(ROIs', means(:,1), means(:,2), means(:,3), 'VariableNames', {'ROIs', 'visual', 'auditory', 'supramodal'});

figure;
bar(means);
xticks(1:N_ROIs);
xticklabels(replace(ROIs, '_', ' '));
ylabel('PSC Ratio (WM/WM+SMC)');
legend({'Visual', 'Auditory', 'Supramodal'});
yline(0.2, '--r');

%% Create LME table

LME_table_full = table();
for ss = 1:N
    subjCode = subjCodes{ss};
    for rr = 1:N_ROIs
        ROI = ROIs{rr};
        for hh = 1:2
            hemi = hemis{hh};
            for cc = 1:3
                contrast = contrasts_simplified{cc};
                PSC_ratio = PSC_ratios(ss,rr,hh,cc);
                coord = coordinates_cut(ss,rr,hh,cc);
                if ~isnan(coord) && ~isnan(PSC_ratio)
                    LME_table_full = [LME_table_full; {PSC_ratio, coord, subjCode, ROI, hemi, contrast}];
                end
            end % contrasts
        end % hemis
    end % ROIs
end %subjs

LME_table_full.Properties.VariableNames = {'PSC_ratio', 'coord', 'subjCode', 'ROI', 'hemi', 'contrast'};
LME_table_full.subjCode = categorical(LME_table_full.subjCode);
LME_table_full.hemi = categorical(LME_table_full.hemi);


%% Fit LMEs
lmes = cell(3,1);
betas = nan(3,1);
CIs = nan(3,2);

for cc = 1:3

    LME_table = LME_table_full(strcmp(LME_table_full.contrast, contrasts_simplified{cc}),:);
    LME_table.coord = (LME_table.coord - coord_means(cc))/coord_stds(cc); % normalize predictor for beta interpretability
    lme = fitglme(LME_table, 'PSC_ratio ~ 1 + coord + (1 + coord | subjCode) + (1 + coord | hemi)')
    lmes{cc} = lme;
    
    residuals = lme.residuals;
    skew = skewness(residuals);
    kurt = kurtosis(residuals);
    disp([contrasts_simplified{cc} ' skewness: ' num2str(skew), '| kurtosis: ' num2str(kurt)]);
    
    figure; 
    qqplot(residuals);
    title(['QQ Plot ' contrasts_simplified{cc}]);

    figure; 
    fitted_values = predict(lme);
    scatter(fitted_values, residuals);
    title(['Heteroscedasticity Check: ' contrasts_simplified{cc}]);
    xlabel('Fitted Values');
    ylabel('Residuals');

    % Breusch-Pagan test for heteroscedasticity
    n  = height(lme.Variables);
    df = lme.NumPredictors;    
    x = lme.Variables(:,lme.PredictorNames);
    aux = fitlm(x,residuals.^2);
    T = aux.Rsquared.Ordinary*n;
    P = 1-chi2cdf(abs(T),df)

    % Calculate robust covariance (sandwich estimator)
    X = lme.designMatrix;
    meat = zeros(2,2);
    for ss = 1:N % loop through subjs random effect
        idx = ismember(LME_table.subjCode, subjCodes{ss});
        Xg = X(idx,:);
        rg = residuals(idx);
        meat = meat + (Xg' * rg) * (Xg' * rg)';
    end
    
    bread = inv(X' * X);
    covB_robust = bread * meat * bread;
    SEs_robust = sqrt(diag(covB_robust));
    
    % Wald test for significance
    W = betas(2)/SEs_robust(2);
    p = 1 - normcdf(W); % 1-tailed test because we have a directional hypothesis

    % Confidence interval
    betas(cc) = lme.Coefficients.Estimate(2);
    CI(cc,:) = [betas(cc)-1.96*SEs_robust(2),betas(cc)+1.96*SEs_robust(2)];

    % Beta represents the change in PSC ratio for a 1 SD increase in
    % coordinate (this is a good way to express the effect size)

end


%% Plot posterior-anterior coord against PSC ratio
coordinates_cut = coordinates(:,:,:,[1,3,5]);
contrasts_simplified = {'visual', 'auditory', 'supramodal'};
for cc = 1:3
    ratios = PSC_ratios(:,:,:,cc);
    coords = coordinates_cut(:,:,:,cc);
    figure;
    scatter(coords(:), ratios(:), 'filled');
    lsline;
    hold on;
    scatter(mean(coords,[1,3], 'omitnan'), mean(ratios,[1,3], 'omitnan'), 50, 'r', 'filled');
    title([ contrasts_simplified{cc} ' | B=' num2str(round(betas(cc), 3)) ' | CI: [' num2str(round(CI(cc,1),3)) '-' num2str(round(CI(cc,2),3)) ']']);
    xlabel('Posterior->Anterior Coordinate');
    ylabel('PSC Ratio');
    legend({'individual', 'least squares line', 'mean per ROI'})
end


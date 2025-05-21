%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to assess the overlap between supramodal,
% visual, and auditory ROIs using dice coefficients.
% Tom Possidente - Feb 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/projects/spacetime_network/functions/'));
ccc;

%% Initialize Key Variables
subjCodes = {'MM', 'PP', 'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS', 'AI'}; % 21 subjs - union of resting state and localizer analyses
N = length(subjCodes);

ROI_dir = '/projectnb/somerslab/tom/projects/spacetime_network/data/ROIs/';

av_ROI_names = {'sPCS', 'iPCS', 'midIFS', 'pVis', 'tgPCS', 'cIFSG', 'CO', 'FO', 'cmSFG', 'pAud'};
N_av = length(av_ROI_names);
sm_ROI_names = {'sm_sPCS', 'sm_iPCS', 'sm_midFSG', 'sm_aINS', 'sm_dACC', 'sm_preSMA'};
N_sm = length(sm_ROI_names);
hemis = {'lh', 'rh'};
N_hemis = length(hemis);

dice_coeffs = nan(N, N_hemis, N_av, N_sm);

%% Loop through subjs and ROIs and assess overlap 
for nn = 1:N
    subjCode = subjCodes{nn};
    for hh = 1:N_hemis
        hemi = hemis{hh};

        % Load annotation files
        [~, labels_av, ctable_av] = read_annotation([ROI_dir hemi '.' subjCode '_ROIs.annot']); 
        [~, labels_sm, ctable_sm] = read_annotation([ROI_dir hemi '.' subjCode '_smROIs.annot']); 

        for av = 1:N_av
            av_ROI_name = av_ROI_names{av};
            av_ROI_code = ctable_av.table(ismember(ctable_av.struct_names, av_ROI_name),5);
            av_ROI_mask = labels_av == av_ROI_code;
            for sm = 1:N_sm
                sm_ROI_name = sm_ROI_names{sm};
                sm_ROI_code = ctable_sm.table(ismember(ctable_sm.struct_names, sm_ROI_name),5);
                sm_ROI_mask = labels_sm == sm_ROI_code;
                dice_coeffs(nn,hh,av,sm) = (2*sum(av_ROI_mask & sm_ROI_mask)) / (sum(av_ROI_mask) + sum(sm_ROI_mask));
            end
        end
    end
end


%% Plot heatmap of dice coef matrix
groupavg_dice= squeeze(mean(dice_coeffs, 1));
groupavg_dice_lh = squeeze(groupavg_dice(1,:,:));
groupavg_dice_rh = squeeze(groupavg_dice(2,:,:));

figure; 
heatmap(replace(sm_ROI_names, '_', ' '), av_ROI_names, groupavg_dice_lh)
title('lh');

figure; 
heatmap(replace(sm_ROI_names, '_', ' '), av_ROI_names, groupavg_dice_rh)
title('rh');

%% Test if hemispheres are significantly different
ps = nan(N_av, N_sm);
for av = 1:N_av
    for sm = 1:N_sm
        [~,ps(av,sm)] = ttest(squeeze(dice_coeffs(:,1,av,sm)), squeeze(dice_coeffs(:,2,av,sm)));
    end
end

%% Combine hemispheres because they are not significantly different after testing
dice_coeffs_nohem = squeeze(mean(dice_coeffs, 2));

%% Test significant difference from 0
ps = nan(N_av, N_sm);
ts = nan(N_av, N_sm);
sds = nan(N_av, N_sm);
for av = 1:N_av
    for sm = 1:N_sm
        [~,ps(av,sm),~,t] = ttest(squeeze(dice_coeffs_nohem(:,av,sm)));
        ts(av,sm) = t.tstat;
        sds(av,sm) = t.sd;
    end
end

groupavg_nohem = squeeze(mean(dice_coeffs_nohem,1));

figure; 
heatmap(replace(sm_ROI_names, '_', ' '), av_ROI_names, groupavg_nohem)


figure; 
heatmap(replace(sm_ROI_names, '_', ' '), av_ROI_names, ps)
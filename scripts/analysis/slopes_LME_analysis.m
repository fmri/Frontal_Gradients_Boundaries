%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to analyze the ROI t-stat slope data
%%% using LMEM to test for differences by hemisphere, contrast, and
%%% modailty
%%%
%%% Tom Possidente - May 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables
ROI_tstat_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/individual_ROI_tstats/';
ROIs = {'preSMA', 'inf_lat_frontal', 'aINS'};
N_ROIs = length(ROIs);

hemis = {'lh', 'rh'};
N_hemis = length(hemis);

modalities = {'visual', 'auditory'};
N_modalities = length(modalities);

contrasts = {'sensory', 'WM'};
N_contrasts = length(contrasts);

subjCodes = {'MK','AB', 'AD','LA','AE','TP','NM','AF','AG','AI','GG','UV','PQ','KQ','LN','RT','PT','PL','NS','SL'};
N_subjs = length(subjCodes);

%% Load and unpack slope data
slope_data = table();

count = 0;
for rr = 1:N_ROIs
    ROI = ROIs{rr};
    for mm = 1:N_modalities
        modality = modalities{mm};
        load([ROI_tstat_dir ROI '_' modality '_slopes.mat'], 'slopes');
        for hh = 1:N_hemis
            hemi = hemis{hh};
            for cc = 1:N_contrasts
                contrast = contrasts{cc};
                for ss = 1:N_subjs
                    count = count + 1;
                    subj = subjCodes{ss};
                    slope = slopes(ss, cc, hh);
                    if ~isnan(slope)
                        slope_data(end+1,:) = {abs(slope), hemi, contrast, modality, subj, ROI};
                    end
                end
            end
        end
    end
end

slope_data.Properties.VariableNames = {'slope', 'hemisphere', 'contrast', 'modality', 'subject', 'ROI'};
slope_data.subject = categorical(slope_data.subject);
slope_data.hemisphere = categorical(slope_data.hemisphere);
slope_data.contrast = categorical(slope_data.contrast);
slope_data.modality = categorical(slope_data.modality);
slope_data.ROI = categorical(slope_data.ROI);

%% Fit LMEM
lme = fitglme(slope_data, 'slope ~ 1 + contrast*modality*hemisphere + (1 + contrast*modality*hemisphere | subject) + (1 + contrast*modality*hemisphere | ROI)',...
    'CovariancePattern', {'Diagonal', 'Diagonal'})

emm = emmeans(lme,'unbalanced');
emm.table

plot_psc_emmeans(sortrows(emm.table,'Row','descend'));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to calculate the mean difference in
%%% T-statistics across each ROI using the fitted gradient models.
%%%
%%% Tom Possidente - September 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables
subjCodes = {'MK','AB', 'AD','LA','AE','TP','NM','AF','AG','AI','GG','UV','PQ','KQ','LN','RT','PT','PL','NS','SL'};
N = length(subjCodes);

modelfits_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/Individual_ROI_models/';

hemis = {'lh', 'rh'};
N_hemis = length(hemis);

modalities = {'auditory', 'visual'};
N_modalities = length(modalities);

ROIs = {'preSMA', 'inf_lat_frontal', 'aINS', 'sup_lat_frontal', 'midIFS'};
N_ROIs = length(ROIs);

T_diff = nan(N,N_hemis,N_modalities,N_ROIs);

%% Loop through subjs/ROIs/modalities

for mm = 1:N_modalities
    modality = modalities{mm};
    for rr = 1:N_ROIs
        ROI = ROIs{rr};
        data = load([modelfits_dir ROI '_' modality '_hinge_boundary_models.mat']);
        for hh = 1:N_hemis
            for ss = 1:N
                subjCode = subjCodes{ss};
                curr_data = data.hinge_means_all{ss,hh};
                if ~isempty(curr_data)
                    T_diff(ss,hh,mm,rr) = curr_data.p4 - curr_data.p3;
                end
            end
        end
    end
end

%% 
T_diff(T_diff<=0) = nan;

squeeze(mean(T_diff,[1,2],'omitnan'))

mean(T_diff, 'all', 'omitnan')

squeeze(std(T_diff,[],[1,2],'omitnan'))

std(T_diff, [], 'all', 'omitnan')
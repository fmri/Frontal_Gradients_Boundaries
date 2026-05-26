%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to analyze potential differences between
%%% visual and auditory axes of greatest change within each ROI to look for
%%% directional differences that can be mapped back into RAS space
%%%
%%% Tom Possidente - May 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables
ROIs = {'preSMA', 'inf_lat_frontal', 'aINS'};
N_ROIs = length(ROIs);
subjCodes = {'MK','AB', 'AD','LA','AE','TP','NM','AF','AG','AI','GG','UV','PQ','KQ','LN','RT','PT','PL','NS','SL'};
N_subjs = length(subjCodes);
contrasts = {'auditory', 'visual', 'supramodal'};
N_contrasts = length(contrasts);
hemis = {'lh', 'rh'};
N_hemis = length(hemis);

label_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/probabilistic_allROIs/';


%% Loop over subjs/ROI/contrasts/hemis and extract RAS axis angle of gradient

for hh = 1:N_hemis
    hemi = hemis{hh};
    for rr = 1:N_ROIs
        ROI = ROIs{rr};
        patch = read_patch([label_dir hemi '.' ROI '_prob_thresh5_flat.patch']);
        for ss = 1:N_subjs
            subjCode = subjCodes{ss};
            for cc = 1:N_contrasts
                contrast = contrasts{cc};
            end
        end
    end
end
% For each ROI we need to find 2 vertices that are along the axis of
% greatest change for that subj/contrast. We can then  use the indices of those vertices 
% (from patch.ind) to look up the vertex in RAS coordinates (from label
% file). 
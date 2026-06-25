%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to take individuals surface labels, convert 
% them to volume, then call fs_labels2annot to convert them into a single 
% annotation file. Then convert the annot files into Conn-compatable nii files.
%
% Tom Possidente June 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'))
ccc;

%% Set key variables
modality = 'visual'; % visual, auditory, supramodal
subjCodes = {'MK','AB', 'AD','LA','AE','TP','NM','AF','AG','AI','GG','UV','KQ','LN', 'PL', 'PT','NS','SL'};
ROIs = {'preSMA', 'inf_lat_frontal', 'aINS', 'midIFS', 'sup_lat_frontal'};
N_ROIs = length(ROIs);
recon_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/';
hemis = {'lh', 'rh'};
ROI_dirs = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/subj_specific_01/';

%% Loop through subjs/hemis and create annot files
for ss = 1:length(subjCodes)
    subjCode = subjCodes{ss};
    annotname = [subjCode '_vol_ROIs'];
    for hh = 1:length(hemis)
        ROI_paths = {};
        hemi = hemis{hh};
        for rr = 1:N_ROIs
            ROI_path = [ROI_dirs hemi '.' ROIs{rr} '_' modality '_' subjCode '.label'];
            unix(['mri_label2vol ' ])
            ROI_paths{end+1} = [recon_dir subjCode ROI_dirs{ROI_ind} hemi ROI_hemi_prefix{ROI_ind} ROIs{rr} '.label'];
        end
        fs_labels2annot(ROI_paths, subjCode, hemi, annotname, './', true, true);

        % Convert annot files to Conn-formatted nii files 
        if hh==2 % this function will find the lh annot and combine into 1 so we only need to run it once per subj
            conn_annot2nii_mod([recon_dir subjCode '/label/' hemi '.' annotname '_fsaverage.annot']);
        end
    end
end

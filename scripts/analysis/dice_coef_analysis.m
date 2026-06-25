%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to calculate dice coefficients between 
%%% sensory and WM activation for each ROI in each subj for auditory and 
%%% contrasts.
%%%
%%% Tom Possidente - June 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables
hemis = {'lh', 'rh'};
num_hemis = length(hemis);

subjCodes = {'MK','AB', 'AD','LA','AE','TP','NM','AF','AG','AI','GG','UV','PQ','KQ','LN','RT','PT','PL','NS','SL'};
N_subjs = length(subjCodes);

contrasts = {'vA-vP', 'vP-f'};
keyword = 'visual';

ROIs = {'aMFG', 'midIFS', 'aINS', 'preSMA', 'inf_lat_frontal', 'sup_lat_frontal', 'aIPS', 'VOT', 'cIPS', 'LOT', 'pIPS', 'VO', 'DO'}; % visual
%ROIs = {'aMFG', 'midIFS', 'aINS', 'preSMA', 'inf_lat_frontal','sup_lat_frontal', 'parietal_opercular', 'aIPS', 'ms_post_STSG', 'cIPS'}; % auditory
%ROIs = {'aINS', 'preSMA', 'inf_lat_frontal', 'sup_lat_frontal', 'aIPS', 'cIPS'}; % supramodal
N_ROIs = length(ROIs);

pval_thresh = -log10(0.05);

ROI_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/subj_specific_01/';
data_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';

dice_coefs = nan(N_subjs, N_ROIs, num_hemis);

%% Loop over subjs/hemis and calculate dice coef between contrasts

for ss = 1:N_subjs
    subjCode = subjCodes{ss};
    for hh = 1:num_hemis
        hemi = hemis{hh};
        for rr = 1:N_ROIs
            ROI_name = ROIs{rr};
            label_file = [ROI_dir hemi '.' ROI_name '_' keyword '_' subjCode '.label'];
            if ~isfile(label_file)
                continue;
            end
            ROI = readtable(label_file, "FileType","text");
            total_sig = 0;
            vertices = true(height(ROI),1);
            for cc = 1:2
                if endsWith(contrasts{cc}, '-f')
                    contrast = ['f-' contrasts{cc}(1:2)];
                    pmap = MRIread([data_dir subjCode '/localizer/localizer_contrasts_sm0_' hemi '/' contrast '/sig.nii.gz']);
                    pmap.vol = -pmap.vol;
                else
                    pmap = MRIread([data_dir subjCode '/localizer/localizer_contrasts_sm0_' hemi '/' contrasts{cc} '/sig.nii.gz']);
                end
                pmap_inROI = pmap.vol(ROI.Var1+1);
                vertices(pmap_inROI<pval_thresh) = false;
                total_sig = total_sig + sum(pmap_inROI>=pval_thresh);
            end
            dice_coefs(ss,rr,hh) = 2*sum(vertices) / total_sig;
        end
    end
end

%% Look at means
means = squeeze(mean(dice_coefs,1,'omitnan'));
sds = squeeze(std(dice_coefs,1,'omitnan'));
means_table = table(ROIs', means(:,1), means(:,2));
means_table.Properties.VariableNames = {'ROIs', 'lh', 'rh'}

figure;
bar(1:N_ROIs, [means_table.lh, means_table.rh]);
xticklabels(ROIs);
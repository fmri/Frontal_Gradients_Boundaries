%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to compute visual and auditory
%%% activation separately for each ROI across the cortex.
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

contrasts = {'aAaP-f', 'vAvP-f'};
contrast_types = {'visaud', 'visaud'};
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
            
            data = MRIread([data_dir subjCode '/localizer/localizer_contrasts_' hemi '/' contrast '/cespct.nii.gz']);

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
PSC_ratios = PSCs(:,:,:,2) ./ (PSCs(:,:,:,1)+PSCs(:,:,:,2)); % vis / vis+aud

means = squeeze(mean(PSC_ratios, [1,3], 'omitnan'));

figure;
bar(means);
xticks(1:N_ROIs);
xticklabels(replace(ROIs, '_', ' '));
ylabel('PSC Ratio (vis/aud+vis)');
yline(0.2, '--r');
yline(0.8, '--r');

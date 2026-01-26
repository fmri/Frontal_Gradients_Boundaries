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
ROIs = {'aINS', 'aIPS', 'aMFG', 'ant_temporal', 'cIPS', 'CO', 'DO', 'inf_aINS', 'iPCS', 'LOT', 'midIFS', 'midINS', ...
        'ms_post_STSG', 'MT', 'parietal_opercular', 'pIPS', 'post_col_sulc', 'post_temporal', 'preSMA', 'sPCS',...
        'sup_aINS', 'tgPCS', 'VO', 'VOT'};
N_ROIs = length(ROIs);

hemis = {'lh', 'rh'};

subjCodes = {'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS', 'AI'};
N = length(subjCodes);

contrasts = {'vP-f', 'vA-vP', 'aP-f', 'aA-aP', 'aPvP-f', 'vAaA-vPaP'};
contrast_types = {'visual', 'visual', 'auditory', 'auditory', 'supramodal', 'supramodal'};
N_contrasts = length(contrasts);

N_vertices = 163842;
vertex_inds = 1:N_vertices;

ROI_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/subj_specific/';
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

%% Plot posterior-anterior coord against PSC ratio
coordinates_cut = coordinates(:,:,:,[1,3,5]);
titles = {'visual', 'auditory', 'supramodal'};
for cc = 1:3
    ratios = PSC_ratios(:,:,:,cc);
    coords = coordinates_cut(:,:,:,cc);
    figure;
    scatter(coords(:), ratios(:), 'filled');
    title([ titles{cc} ' r=' num2str(round(corr(ratios(:), coords(:),'rows','complete'), 3)) ]);
end

%%
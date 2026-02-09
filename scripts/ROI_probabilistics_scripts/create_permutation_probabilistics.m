%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to create probabilistic maps from
%%% randomly permuted groups of WM and SMC contrasts in order to create a null
%%% distribution for eventual hypothesis testing with true WM and SMC
%%% contrast groups
%%%
%%% Tom Possidente - January 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/');
ccc;

tic; 

%% Initialize Key Variables

iterations = 10000; 
contrasts = {'vA-vP', 'aA-aP'; 'f-vP', 'f-aP'};
keyword= 'supramodal'; % used for naming nii files (auditory or visual)
contrast_labels = {'WM', 'SMC'};
N_contrasts = length(contrasts);

subjCodes = {'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS', 'AI'};
N = length(subjCodes);

N_vertices = 163842;
hemis = {'lh', 'rh'};
N_hemis = length(hemis);
pval_thresh = -log10(0.05); % pvals in freesurfer files are expressed as -log10(p)

permutations = rand(N, iterations)>0.5; % swap or not for each subj in each iteration (50% chance)

data_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';
permuted_maps_outdir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/permutation_probabilistic_maps_WMSMC/';

%% Create alternative contrast names if needed
 flip = zeros(N_contrasts);
for cc = 1:N_contrasts
    contrast = contrasts{cc};
    if length(contrasts{cc})==4 & strcmp(contrast(end-2:end), 'P-f')
        split_str = split(contrast,'-');
        contrasts_alt{cc} = [split_str{2} '-' split_str{1}]; % contrast string coded backwards for single modality SMC contrasts
        flip(cc) = 1;
    else
        contrasts_alt{cc} = contrast;
    end
end

%% Load data from each contrast for each subj
stat_data = nan(N_vertices, N_hemis, N_contrasts, N);

for hh = 1:N_hemis
    hemi = hemis{hh};
    for nn = 1:N
        subjcode = subjCodes{nn};
        if strcmp(keyword, 'supramodal') 
            for cc = 1:N_contrasts
                data_temp = nan(N_vertices, 2);
                for ss = 1:2
                    contrast = contrasts{cc,ss};
                    data = MRIread([data_dir subjcode '/localizer/localizer_contrasts_' hemi '/' contrast '/sig.nii.gz']);
                    if length(contrast)==4 & strcmp(contrast(1:2), 'f-')
                        data.vol = -data.vol; % if contrast was backwards, flip it to be correct
                    end
                    data_temp(:,ss) = data.vol>pval_thresh;
                end
                stat_data(:,hh,cc,nn) = data_temp(:,1) & data_temp(:,2);
            end
        else
            for cc = 1:N_contrasts
                contrast = contrasts{cc};
                data = MRIread([data_dir subjcode '/localizer/localizer_contrasts_' hemi '/' contrasts_alt{cc} '/sig.nii.gz']);
                if flip(cc)
                    data.vol = -data.vol; % if contrast was backwards, flip it to be correct
                end
                stat_data(:,hh,cc,nn) = data.vol>pval_thresh;
            end
        end
    end
end

%% Create probabilistic map for each permutation WM and SMC labels 
grand_probmap = zeros(N_vertices, N_hemis);

for ii = 1:iterations

    % Get permuted data in 2 groups
    current_perm = permutations(:,ii); % get N preallocated random 1s and 0s

    subj_data1 = squeeze(stat_data(:,:,1,:)); % select all WM data
    subj_data1(:,:,current_perm==1) = squeeze(stat_data(:,:,2,current_perm==1)); % replace WM with SMC data where perm is 1

    subj_data2 = squeeze(stat_data(:,:,2,:)); % select all SMC data
    subj_data2(:,:,current_perm==1) = squeeze(stat_data(:,:,1,current_perm==1)); % replace SMC with WM data where perm is 1

    % Create probabilistic for each group
    probmap1 = sum(subj_data1,3); % summing over subjs 
    probmap2 = sum(subj_data2,3);

    % Add to grand_probmap
    grand_probmap = grand_probmap + probmap1 + probmap2; 

    % Save out lh probabilistics as nii files
    data.vol = probmap1(:,1)'; % reusing "data" structure and overwriting vol data so that all the other fields are prefilled and correct
    MRIwrite(data, [permuted_maps_outdir, 'lh_iter' num2str(ii) '_group1_' keyword '.nii']);
    data.vol = probmap2(:,1)';
    MRIwrite(data, [permuted_maps_outdir, 'lh_iter' num2str(ii) '_group2_' keyword '.nii']);

    % Save out rh probabilistic as nii files
    data.vol = probmap1(:,2)'; % reusing "data" structure and overwriting vol data so that all the other fields are prefilled and correct
    MRIwrite(data, [permuted_maps_outdir, 'rh_iter' num2str(ii) '_group1_' keyword '.nii']);
    data.vol = probmap2(:,2)';
    MRIwrite(data, [permuted_maps_outdir, 'rh_iter' num2str(ii) '_group2_' keyword '.nii']);

    if mod(ii,500)==0
        disp(num2str(round(ii/iterations,2)));
    end

end

% Save out grand probmap
data.vol = grand_probmap(:,1)';
MRIwrite(data, [permuted_maps_outdir, 'lh_grand_probmap_' keyword '.nii']);
data.vol = grand_probmap(:,2)';
MRIwrite(data, [permuted_maps_outdir, 'rh_grand_probmap_' keyword '.nii']);
grand_probmap_thresh = iterations*N_contrasts*N * 0.1;

%% Create probabilistic map for unpermuted data

WM_data = squeeze(stat_data(:,:,1,:));
SMC_data = squeeze(stat_data(:,:,2,:));

WM_probmap = sum(WM_data, 3);
SMC_probmap = sum(SMC_data, 3);

data.vol = WM_probmap(:,1)'; % reusing "data" structure and overwriting vol data so that all the other fields are prefilled and correct
MRIwrite(data, [permuted_maps_outdir, 'lh_original_WM_' keyword '.nii']);
data.vol = SMC_probmap(:,1)';
MRIwrite(data, [permuted_maps_outdir, 'lh_original_SMC_' keyword '.nii']);

data.vol = WM_probmap(:,2)'; % reusing "data" structure and overwriting vol data so that all the other fields are prefilled and correct
MRIwrite(data, [permuted_maps_outdir, 'rh_original_WM_' keyword '.nii']);
data.vol = SMC_probmap(:,2)';
MRIwrite(data, [permuted_maps_outdir, 'rh_original_SMC_' keyword '.nii']);



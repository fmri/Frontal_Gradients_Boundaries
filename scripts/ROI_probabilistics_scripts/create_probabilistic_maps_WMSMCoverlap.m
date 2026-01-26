%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to create a probabilistic map of the
%%% overlap between WM and SMC contrasts
%%% Tom Possidente - December 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/');
ccc;

%% Get Subj Codes

subjCodes = {'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS', 'AI'};
N = length(subjCodes);


%% Initialize variables
contrasts = {'aA-aP', 'f-aP'};
keyword = 'aud';

data_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';
outdir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/';

p_thresh = 1.3;
N_contrasts = length(contrasts);
hemis = {'lh', 'rh'};
N_hemis = length(hemis);
N_vertices = 163842;

%% Loop through subjs, extract contrast stat data, make probabilistic nii
for hh = 1:N_hemis
    hemi = hemis{hh};
    subj_vertices.vol = zeros(1,N_vertices);
    for ss = 1:N
        subjCode = subjCodes{ss};
        binarized_stat_data = true(1, N_vertices);
        for cc = 1:N_contrasts

            stat_path = [data_dir subjCode '/localizer/localizer_contrasts_' hemi '/' contrasts{cc} '/sig.nii.gz'];
            stat_data = MRIread(stat_path);

            if ismember(contrasts{cc}, {'f-vP', 'f-aP', 'f-tP'})
                stat_data.vol = -stat_data.vol; % reverse contrast for interpretability
            end

            binarized_stat_data = stat_data.vol >= p_thresh & binarized_stat_data;
            
        end

        subj_vertices.vol = subj_vertices.vol + binarized_stat_data;

    end

    % Save nii file
    MRIwrite(subj_vertices, [outdir hemi '_' keyword 'WMSMC_overlap_probabilistic.nii']);
    disp(['Finishing ' hemi ' ' contrasts]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to create a probabilistic map of 
%%% the individual subject difference maps between the WM and SMC contrasts 
%%% Tom Possidente - April 206
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/');
ccc;

%% Get Subj Codes

subjCodes = {'MK', 'AB', 'AD', 'LA', 'AE', 'TP', 'NM', 'AF', 'AG', 'GG', 'UV', 'PQ', 'KQ', 'LN', 'RT', 'PT', 'PL', 'NS', 'AI'};
N = length(subjCodes);


%% Initialize variables
contrasts = {'vA-vP', 'f-vP'};
keyword = 'vis';

data_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';
outdir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/probabilisitc_maps/';

N_contrasts = length(contrasts);
hemis = {'lh', 'rh'};
N_hemis = length(hemis);
N_vertices = 163842;

%% Loop through subjs, extract contrast stat data, make probabilistic nii

for hh = 1:N_hemis
    hemi = hemis{hh};
    probabilistic = zeros(1,N_vertices);
    for ss = 1:N
        subjCode = subjCodes{ss};

        stat_path1 = [data_dir subjCode '/localizer/localizer_contrasts_' hemi '/' contrasts{1} '/t.nii.gz'];
        stat_data1 = MRIread(stat_path1);
        if ismember(contrasts{1}, {'f-vP', 'f-aP', 'f-tP'})
            stat_data1.vol = -stat_data1.vol; % reverse contrast for interpretability
        end
        
        stat_path2 = [data_dir subjCode '/localizer/localizer_contrasts_' hemi '/' contrasts{2} '/t.nii.gz'];
        stat_data2 = MRIread(stat_path2);
        if ismember(contrasts{2}, {'f-vP', 'f-aP', 'f-tP'})
            stat_data2.vol = -stat_data2.vol; % reverse contrast for interpretability
        end

        stat_data1.vol(stat_data1.vol<0) = 0;
        stat_data2.vol(stat_data2.vol<0) = 0;
        stat_diffs = stat_data1.vol - stat_data2.vol;
        t_thresh = prctile(stat_diffs, [10,90])
        stat_diffs_WM = stat_diffs>t_thresh(2);
        stat_diffs_SMC = stat_diffs<t_thresh(1);
        probabilistic = probabilistic + stat_diffs_WM;
        probabilistic = probabilistic - stat_diffs_SMC;
    end

    % Save nii file
    stat_data1.vol = probabilistic;
    MRIwrite(stat_data1, [outdir hemi '_' keyword '_WMSMC_diff_probabilistict.nii']);
    disp(['Finishing ' hemi]);

end

%% Supramodal case
contrasts1 = {'vA-vP', 'aA-aP'};
contrasts2 = {'f-vP', 'f-aP'};
keyword = 'supra';

for hh = 1:N_hemis
    hemi = hemis{hh};
    probabilistic = zeros(1,N_vertices);
    for ss = 1:N
        subjCode = subjCodes{ss};
        for cc = 1:2
            stat_path1 = [data_dir subjCode '/localizer/localizer_contrasts_' hemi '/' contrasts1{cc} '/t.nii.gz'];
            stat_data1 = MRIread(stat_path1);
            if ismember(contrasts1{cc}, {'f-vP', 'f-aP', 'f-tP'})
                stat_data1.vol = -stat_data1.vol; % reverse contrast for interpretability
            end
            
            stat_path2 = [data_dir subjCode '/localizer/localizer_contrasts_' hemi '/' contrasts2{cc} '/t.nii.gz'];
            stat_data2 = MRIread(stat_path2);
            if ismember(contrasts2{cc}, {'f-vP', 'f-aP', 'f-tP'})
                stat_data2.vol = -stat_data2.vol; % reverse contrast for interpretability
            end            
            stat_data1.vol(stat_data1.vol<0) = 0;
            stat_data2.vol(stat_data2.vol<0) = 0;
            stat_diffs = stat_data1.vol - stat_data2.vol;
            t_thresh = prctile(stat_diffs, [10,90])
            stat_diffs_WM{cc} = stat_diffs>t_thresh(2);
            stat_diffs_SMC{cc} = stat_diffs<t_thresh(1);
        end

        probabilistic = probabilistic + (stat_diffs_WM{1} & stat_diffs_WM{2});
        probabilistic = probabilistic - (stat_diffs_SMC{1} & stat_diffs_SMC{2});
    end

    % Save nii file
    stat_data1.vol = probabilistic;
    MRIwrite(stat_data1, [outdir hemi '_' keyword '_WMSMC_diff_probabilistict.nii']);
    disp(['Finishing ' hemi]);

end


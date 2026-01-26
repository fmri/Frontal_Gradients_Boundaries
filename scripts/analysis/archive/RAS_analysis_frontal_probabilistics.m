%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to calculate the center of mass of
%%% probabilistic ROIs for WM and SMC contrasts, then compare those RAS
%%% locations (direction and distance)
%%%
%%% Tom Possidente - November 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables

ROI_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/';

hemis = {'lh', 'rh'};
ROI_names = {'sPCS', 'iPCS', 'midIFS', 'aINS', 'preSMA'};
contrasts = {'vA-vP', 'vP-f', 'aA-aP', 'aP-f'};

N_hemis = length(hemis);
N_ROIs = length(ROI_names);
N_contrasts = length(contrasts);

COMs = nan(N_ROIs, N_contrasts, N_hemis, 3); % storage for center of masses RAS coordinates

%% Loop through ROI probabilistic data and calculate center of mass for each

for hh = 1:N_hemis
    hemi = hemis{hh};
    for cc = 1:N_contrasts
        contrast = contrasts{cc};

        if strcmp(contrast(end-1:end), '-f')
            contrast_str_mod = strsplit(contrast, '-');
            contrast_str_nii = [contrast_str_mod{2} '-' contrast_str_mod{1}];
        else
            contrast_str_nii = contrast;
        end

        % Load probabilistic data for full hemisphere
        prob_data_path = [ROI_dir hemi '_' contrast_str_nii '_intersection_probabilistic_visualization.nii'];
        prob_data = MRIread(prob_data_path);

        for rr = 1:N_ROIs
            ROI_name = ROI_names{rr};

            % Load label file
            label_path = [ROI_dir hemi '_' ROI_name '_probabilistic_' replace(contrast, '-', '_') '_thresh5.label'];
            label = readtable(label_path, 'FileType','text');

            % Get #subjs per vertex in ROI
            ROI_data = prob_data.vol(label{:,1}+1); 

            % Calculate center of mass using RAS coords and #subjs per vertex 
            COMs(rr,cc,hh,:) = ROI_data * label{:,2:4} / sum(ROI_data);

        end
    end
end

%% Compare direction and distance between COMs
compare = [1,2; 3,4];
compare_names = {'vis_WMvsSMC', 'aud_WMvsSMC'};
coord_diffs_RAS = table();

counter = 0;
for hh = 1:N_hemis
    for cc = 1:size(compare,1)
        for rr = 1:N_ROIs
           counter = counter + 1; 
           COM_WM = COMs(rr,compare(cc,1),hh,:);
           COM_SMC = COMs(rr,compare(cc,2),hh,:);
           % figure;
           % scatter3(COM_WM(1), COM_WM(2), COM_WM(3)); hold; 
           % scatter3(COM_SMC(1), COM_SMC(2), COM_SMC(3));
           % xlabel('R'); ylabel('A'); zlabel('S');
           % title([hemis{hh} ' ' ROI_names{rr}]);
           % legend({'WM', 'SMC'});

            coord_diffs_RAS{counter,1} = hemis{hh};
            coord_diffs_RAS{counter,2} = compare_names{cc};
            coord_diffs_RAS{counter,3} = ROI_names(rr);
            coord_diffs_RAS{counter,4} = squeeze(COM_WM-COM_SMC)';
            coord_diffs_RAS{counter,5} = pdist(squeeze([COM_WM; COM_SMC]));
        end
    end
end

coord_diffs_RAS.Properties.VariableNames = {'hemisphere', 'contrast_comparison', 'ROI', 'RAS_difference', 'distance(mm)'}; 



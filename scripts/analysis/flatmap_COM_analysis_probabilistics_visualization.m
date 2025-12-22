%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to calculate the center of mass of
%%% probabilistic ROIs for WM and SMC contrasts, then compare those flatmap
%%% XY locations (direction and distance)
%%%
%%% Tom Possidente - November 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables

ROI_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/';

hemis = {'lh', 'rh'};
ROI_names = {'sPCS', 'iPCS', 'midIFS', 'aINS', 'preSMA'};
% contrasts = {'vA-vP', 'vP-f', 'aA-aP', 'aP-f', 'vP-f_ap-f', 'vA-vP_aA-aP'};
% set_names = {'visWMSMC', 'visWMSMC', 'audWMSMC', 'audWMSMC', 'WMSMC', 'WMSMC'};
contrasts = {'vA-vP', 'vP-f'};
set_names = {'visWMSMC', 'visWMSMC'}; % WMSMC, visWMSMC. audSMC

N_hemis = length(hemis);
N_ROIs = length(ROI_names);
N_contrasts = length(contrasts);

COMs = nan(N_ROIs, N_contrasts, N_hemis, 2); % storage for center of masses XY coordinates
COM_vert_inds = nan(N_ROIs, N_contrasts, N_hemis); % storage for vertex index closest to COM

plotting_diagnostics = true;

%% Loop through ROI probabilistic data and calculate center of mass for each

for hh = 1:N_hemis
    hemi = hemis{hh};
    for cc = 1:N_contrasts
        contrast = contrasts{cc};
        set_name = set_names{cc};

        % Load probabilistic data for full hemisphere
        prob_data_path = [ROI_dir hemi '_' contrast '_intersection_probabilistic_visualization.nii'];
        prob_data = MRIread(prob_data_path);

        for rr = 1:N_ROIs
            ROI_name = ROI_names{rr};

            % Load patch file
            patch_path = [ROI_dir hemi '.' ROI_name '_probabilistic_' set_name '_thresh5_flat.patch'];
            if strcmp(ROI_name, 'aINS') & strcmp(hemi, 'rh') & strcmp(set_name, 'WMSMC') % missing aINS for SMC condition
                COM_vert_inds(rr,cc,hh) = nan;
                COMs(rr,cc,hh,:) = nan;
                continue
            else
                patch = read_patch(patch_path);
            end

            % Load label file
            label_path = [ROI_dir hemi '_' ROI_name '_probabilistic_' replace(contrast, '-', '_') '_thresh5.label'];
            label = readtable(label_path, 'FileType', 'text');

            % Get #subjs per vertex in ROI
            inds_label = label{:,1};
            patch_in_label = ismember(inds_label, patch.ind);
            perc_labelinpatch = sum(patch_in_label)/length(inds_label);
            disp(perc_labelinpatch);
            inds = inds_label(patch_in_label);
            inds = sort(inds);
            ROI_data = prob_data.vol(inds+1);

            % Calculate center of mass using XY coords and #subjs per vertex
            ind_mask = ismember(patch.ind, inds);
            x = patch.x(ind_mask);
            y = patch.y(ind_mask);
            patch_inds = patch.ind(ind_mask);
            COMs(rr,cc,hh,:) = ROI_data * [x; y]' / sum(ROI_data);

            % Find and record vertex index closest to COM point
            [vert_dist, closest_ind] = min( sqrt( sum( ([x;y] - squeeze(COMs(rr,cc,hh,:))).^2 ) ) ); 
            if vert_dist>0.999
                keyboard;
            end

            COM_vert_inds(rr,cc,hh) = patch_inds(closest_ind);
            
            if plotting_diagnostics
                figure; 
                scatter(x,y,[],ROI_data, 'filled'); hold on; 
                scatter(COMs(rr,cc,hh,1), COMs(rr,cc,hh,2), 100, 'r', 'filled'); 
                scatter(x(closest_ind), y(closest_ind), 100, 'o', 'filled');
                title([hemi ' ' ROI_name ' ' contrast])
            end

        end
    end
end

%% Compare direction and distance between COMs
% compare = [1,2; 3,4];
% compare_names = {'vis_WMvsSMC', 'aud_WMvsSMC'};
compare = [1,2];
compare_names = {'WMvsSMC'};
coord_diffs = table();

counter = 0;
for hh = 1:N_hemis
    hemi = hemis{hh};
    for cc = 1:size(compare,1)
        for rr = 1:N_ROIs
            ROI_name = ROI_names{rr};
            counter = counter + 1;
            COM_WM = squeeze(COMs(rr,compare(cc,1),hh,:));
            COM_SMC = squeeze(COMs(rr,compare(cc,2),hh,:));
            
            if plotting_diagnostics
                if strcmp(ROI_name, 'aINS') & strcmp(hemi, 'rh') & strcmp(set_name, 'WMSMC') % missing aINS for SMC condition
                    disp('rh aINS unavailable');
                else
                    patch_path = [ROI_dir hemi '.' ROI_name '_probabilistic_' set_names{cc+1} '_thresh5_flat.patch'];
                    patch = read_patch(patch_path);
                    figure; 
                    scatter(patch.x, patch.y); hold on; 
                    scatter(COM_WM(1), COM_WM(2), 25, 'filled');
                    scatter(COM_SMC(1), COM_SMC(2), 25, 'filled');
                    legend({'coordinates', 'WM peak', 'SMC peak'});
                    title([hemis{hh} ' ' ROI_name ' ' set_names{cc+1}]);
                end
            end

            coord_diffs{counter,1} = hemis{hh};
            coord_diffs{counter,2} = compare_names{cc};
            coord_diffs{counter,3} = {ROI_name};
            coord_diffs{counter,4} = squeeze(COM_WM-COM_SMC)';
            coord_diffs{counter,5} = pdist([COM_WM'; COM_SMC']);
        end
    end
end

coord_diffs.Properties.VariableNames = {'hemisphere', 'contrast_comparison', 'ROI', 'XY_difference', 'distance(mm)'};


%% Make labels out of COM points for visualization
label_cortex_lh = readtable([ROI_dir 'lh_cortex.label'], 'FileType','text');
label_cortex_rh = readtable([ROI_dir 'rh_cortex.label'], 'FileType','text');
label_cortex = {label_cortex_lh, label_cortex_rh};

for hh = 1:N_hemis
    hemi = hemis{hh};
    for cc = 1:size(compare,1)
        for rr = 1:N_ROIs

            ROI_name = ROI_names{rr};
            vert_WM = COM_vert_inds(rr,compare(cc,1),hh);
            vert_SMC = COM_vert_inds(rr,compare(cc,2), hh);
                
            new_label = label_cortex{hh}(ismember(label_cortex{hh}{:,1}, [vert_WM, vert_SMC]),:);

            % Make label file
            label_fname = ['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/' hemi '.peaks_' ROI_name '_' compare_names{cc} '.label'];
            label_file = fopen(label_fname,'w');
            fprintf(label_file, ['#!ascii label  , from subject  vox2ras=TkReg\n' num2str(size(new_label,1)) '\n']);
            writetable(new_label, label_fname, 'Delimiter', 'tab', 'WriteMode', 'append', 'FileType', 'text');
            fclose(label_file);

        end
    end
end






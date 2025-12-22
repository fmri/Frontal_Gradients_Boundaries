%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to convert label files into patches and
%%% flatten them using freesurfer's label2patch and mris_flatten functions
%%%
%%% Tom Possidente - December 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize key variables
ROI_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/';
patch_tempdir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/fsaverage/surf/';
hemis = {'lh', 'rh'};
ROI_names = {'sPCS', 'iPCS', 'midIFS', 'aINS', 'preSMA'};


N_hemis = length(hemis);
N_ROIs = length(ROI_names);
N_fnames = N_ROIs*N_hemis;

fnames = {};
for hh = 1:N_hemis
    hemi = hemis{hh};
    for rr = 1:N_ROIs
        ROI_name = ROI_names{rr};
        fnames{end+1} = [hemi '.' ROI_name '_probabilistic_'];
    end
end

%% Loop through hemispheres, ROIs and make patch then flatten it
contrasts = {'vA_vP', 'vP_f'; 'aA_aP', 'aP_f'};
cnames = {'visWMSMC', 'audWMSMC'};

for nn = 1:N_fnames
    hemi = fnames{nn}(1:2);

    for cc = 1:2
        contrast_set = contrasts(cc,:);

        % Combine labels
        label_path1 = [ROI_dir hemi '_' fnames{nn}(4:end) contrast_set{1} '_thresh5.label'];
        label_path2 = [ROI_dir hemi '_' fnames{nn}(4:end) contrast_set{2} '_thresh5.label'];
        label_merged = [ROI_dir fnames{nn} cnames{cc} '.label'];
        if ~isfile(label_merged)
            unix(['mri_mergelabels -i ' label_path1 ' -i ' label_path2 ' -o ' label_merged])
        end

        % convert label to patch
        patch_outpath = [patch_tempdir fnames{nn} cnames{cc} '_thresh5.patch'];
        if ~isfile(patch_outpath)
            unix(['label2patch -surf white fsaverage ' hemi ' ' label_merged ' ' patch_outpath]);
        end

        % convert patch to flattened patch
        flatpatch_outpath = [ROI_dir fnames{nn} cnames{cc} '_thresh5_flat.patch'];
        if ~isfile(flatpatch_outpath)
            unix(['mris_flatten ' patch_outpath ' ' flatpatch_outpath]);
        else
            continue;
        end

        % Remove intermediate files created during patch converstion/flattening
        unix(['rm ' patch_outpath]);
        unix(['rm ./' fnames{nn} cnames{cc} '_thresh5_flat.patch.out'])
    end
end



%% Do the same for the supramodal versions
contrasts = {'vP_f_aP_f', 'vA_vP_aA_aP'};
for nn = 1:N_fnames
    hemi = fnames{nn}(1:2);

    % Combine labels
    label_path1 = [ROI_dir hemi '_' fnames{nn}(4:end) contrasts{1} '_thresh5.label'];
    label_path2 = [ROI_dir hemi '_' fnames{nn}(4:end) contrasts{2} '_thresh5.label'];
    label_merged = [ROI_dir fnames{nn} 'WMSMC.label'];
    if ~isfile(label_merged)
        unix(['mri_mergelabels -i ' label_path1 ' -i ' label_path2 ' -o ' label_merged])
    end

    % convert label to patch
    patch_outpath = [patch_tempdir fnames{nn} 'WMSMC_thresh5.patch'];
    flatpatch_outpath = [ROI_dir fnames{nn} 'WMSMC_thresh5_flat.patch'];
    if ~isfile(flatpatch_outpath)
        unix(['label2patch -surf white fsaverage ' hemi ' ' label_merged ' ' patch_outpath]);
    end

    % convert patch to flattened patch
    if ~isfile(flatpatch_outpath)
        unix(['mris_flatten ' patch_outpath ' ' flatpatch_outpath]);
    else
        continue;
    end

    % Remove intermediate files created during patch converstion/flattening
    unix(['rm ' patch_outpath]);
    unix(['rm ./' fnames{nn} 'WMSMC_thresh5_flat.patch.out'])
end

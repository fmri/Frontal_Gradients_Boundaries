%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to convert label files into patches and
%%% flatten them using freesurfer's label2patch and mris_flatten functions
%%%
%%% Tom Possidente - December 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize key variables
label_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/permutation_analysis_labels_etc/';
patch_tempdir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/fsaverage/surf/';
hemis = {'lh', 'rh'};
ROI_names = {'sPCS', 'iPCS', 'midIFS', 'aINS', 'preSMA'};
keyword = 'supramodal';

N_hemis = length(hemis);
N_ROIs = length(ROI_names);
N_fnames = N_ROIs*N_hemis;

fnames = {};
for hh = 1:N_hemis
    hemi = hemis{hh};
    for rr = 1:N_ROIs
        ROI_name = ROI_names{rr};
        fnames{end+1} = [hemi '.' keyword '_' ROI_name '_grand_probabilistic'];
    end
end

%% Loop through hemispheres, ROIs and make patch then flatten it

for nn = 1:N_fnames
    hemi = fnames{nn}(1:2);

    % convert label to patch
    label_path = [label_dir fnames{nn} '.label'];
    patch_outpath = [patch_tempdir fnames{nn} '_thresh5.patch'];
    if ~isfile(patch_outpath)
        unix(['label2patch -surf white fsaverage ' hemi ' ' label_path ' ' patch_outpath]);
    end

    % convert patch to flattened patch
    flatpatch_outpath = [label_dir fnames{nn} '_thresh5_flat.patch'];
    if ~isfile(flatpatch_outpath)
        unix(['mris_flatten ' patch_outpath ' ' flatpatch_outpath]);
    else
        continue;
    end

    % Remove intermediate files created during patch converstion/flattening
    unix(['rm ' patch_outpath]);
    unix(['rm ./' fnames{nn} '_thresh5_flat.patch.out'])
end




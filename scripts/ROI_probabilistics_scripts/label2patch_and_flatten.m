%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to convert label files into patches and
%%% flatten them using freesurfer's label2patch and mris_flatten functions
%%%
%%% Tom Possidente - January 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize key variables
ROI_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/probabilistic_allROIs/';
patch_tempdir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/fsaverage/surf/';
hemis = {'lh', 'rh'};
% ROI_names = {'aMFG', 'midIFS', 'aINS', 'preSMA', 'inf_lat_frontal', 'sup_lat_frontal', ...
%              'parietal_opercular', 'aIPS', 'ms_post_STSG', 'VOT', 'cIPS', 'LOT', 'pIPS', 'VO', 'DO'};
ROI_names = {'MT'};

N_hemis = length(hemis);
N_ROIs = length(ROI_names);
N_fnames = N_ROIs*N_hemis;

fnames = {};
for hh = 1:N_hemis
    hemi = hemis{hh};
    for rr = 1:N_ROIs
        ROI_name = ROI_names{rr};
        fnames{end+1} = [hemi '.' ROI_name '_prob_thresh5'];
    end
end

%% Loop through hemispheres, ROIs and make patch then flatten it

for nn = 1:N_fnames
    hemi = fnames{nn}(1:2);

    % convert label to patch
    label_path = [ROI_dir fnames{nn} '.label'];
    patch_outpath = [patch_tempdir fnames{nn} '.patch'];
    if ~isfile(patch_outpath)
        unix(['label2patch -surf white fsaverage ' hemi ' ' label_path ' ' patch_outpath]);
    end

    % convert patch to flattened patch
    flatpatch_outpath = [ROI_dir fnames{nn} '_flat.patch'];
    if ~isfile(flatpatch_outpath)
        unix(['mris_flatten ' patch_outpath ' ' flatpatch_outpath]);
    else
        continue;
    end

    % Remove intermediate files created during patch converstion/flattening
    unix(['rm ' patch_outpath]);
    unix(['rm ./' fnames{nn} '_flat.patch.out'])
end



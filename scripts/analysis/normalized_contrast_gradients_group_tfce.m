%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to produce normalized contrast comparison
% maps comparing WM and sensory activation for 3 modalities
% Tom Possidente - May 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables
fs_nverts = 163842;
stat_pathbase = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/grouplevel/';
tfce_path = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/group_contrast_TFCE/';

hemis = {'lh', 'rh'};
latmed = 'lateral';

frontal_label_lh = ['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/lh.frontal_' latmed '_cortex.label'];
frontal_verts_lh = readtable(frontal_label_lh, 'FileType','text');
frontal_label_rh = ['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/rh.frontal_' latmed '_cortex.label'];
frontal_verts_rh = readtable(frontal_label_rh, 'FileType','text');
frontal_verts = {frontal_verts_lh, frontal_verts_rh};


%% 

contrasts = {'vA-vP', 'f-vP'};
comp_maps = nan(fs_nverts,2);

for hh = 1:2
    hemi = hemis{hh};
    stat_maps = nan(fs_nverts, 2);

    % Get union of significant tfce clusters for each contrast
    contrast1_tfce_p = MRIread([tfce_path contrasts{1} '/' hemi '/_tfce_tstat_fwep_nifti1.nii']);
    contrast1_tfce_sig = contrast1_tfce_p.vol > 1.3;
    contrast2_tfce_p = MRIread([tfce_path contrasts{2} '/' hemi '/_tfce_tstat_fwep_nifti1.nii']);
    contrast2_tfce_sig = contrast2_tfce_p.vol > 1.3;
    union_sig = contrast1_tfce_sig | contrast2_tfce_sig;

    for cc = 1:2 % contrast
        contrast = contrasts{cc};

        % Load tstats
        stat_path = [stat_pathbase hemi '.ces.localizer_groupavg_' contrast '.glmres/osgm/z.mgh'];
        stat_data = MRIread(stat_path);

        if contains(contrast, 'f-')
            stat_data.vol = -stat_data.vol; % reverse contrast for to get sensory drive (passive-fixation)
        end

        % nan out all stats except those in the frontal lobe label
        inlabel_data = stat_data.vol(frontal_verts{hh}.Var1+1);
        stat_data.vol(:) = nan;
        stat_data.vol(frontal_verts{hh}.Var1+1) = inlabel_data;

        % Normalize stats
        disp([hemi ' ' contrast ' stat mean=' num2str(nanmean(stat_data.vol)) ' max=' num2str(max(stat_data.vol)) ])
        neg_tstat_mask = stat_data.vol < 0;
        stat_data.vol(frontal_verts{hh}.Var1+1) = zscore(stat_data.vol(frontal_verts{hh}.Var1+1));
        stat_data.vol(neg_tstat_mask) = nan;
        stat_maps(:,cc) = stat_data.vol;

        % nan out all vertices outside union of significant clusters 
        stat_maps(~union_sig,cc) = nan;

    end

    % Subtract normalized stat maps
    comp_maps(:,hh) = stat_maps(:,1) - stat_maps(:,2); 

    % Save as nii
    outfile = ['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/gradient_niftis/' hemi '_' latmed '_' contrasts{1} '-' contrasts{2} '_groupcomp_frontal_tfce.nii.gz'];
    stat_data.vol = comp_maps(:,hh)';
    MRIwrite(stat_data, outfile)
end


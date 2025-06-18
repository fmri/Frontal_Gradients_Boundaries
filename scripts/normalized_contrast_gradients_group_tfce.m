%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to produce normalized contrast comparison
% maps comparing WM and sensory activation for 3 modalities
% Tom Possidente - May 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables
fs_nverts = 163842;
gstat_pathbase = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/grouplevel/';
tfce_path = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/group_contrast_TFCE/';

hemis = {'lh', 'rh'};
modality = 'v'; % visual = v, auditory = a, tactile = t
latmed = 'medial';

frontal_label_lh = ['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/lh.frontal_' latmed '_cortex.label'];
frontal_verts_lh = readtable(frontal_label_lh, 'FileType','text');
frontal_label_rh = ['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/rh.frontal_' latmed '_cortex.label'];
frontal_verts_rh = readtable(frontal_label_rh, 'FileType','text');
frontal_verts = {frontal_verts_lh, frontal_verts_rh};

%% Loop through subjs, get t-stats for each contrast, normalize them, then subtract them and save as nii

contrasts = {['f-' modality 'P'], [modality 'A-' modality 'P']};
comp_maps = nan(fs_nverts,2);

for hh = 1:2
    hemi = hemis{hh};
    gstat_maps = nan(fs_nverts, 2);

    % Get union of significant tfce clusters for each contrast
    contrast1_tfce_p = MRIread([tfce_path contrasts{1} '/' hemi '/_tfce_tstat_fwep_nifti1.nii']);
    contrast1_tfce_sig = contrast1_tfce_p.vol > 1.3;
    contrast2_tfce_p = MRIread([tfce_path contrasts{2} '/' hemi '/_tfce_tstat_fwep_nifti1.nii']);
    contrast2_tfce_sig = contrast2_tfce_p.vol > 1.3;
    union_sig = contrast1_tfce_sig | contrast2_tfce_sig;

    for cc = 1:2 % contrast
        contrast = contrasts{cc};

        % Load tstats
        gstat_path = [gstat_pathbase hemi '.ces.localizer_groupavg_' contrast '.glmres/osgm/gamma.mgh'];
        gstat_data = MRIread(gstat_path);

        if contains(contrast, 'f-')
            gstat_data.vol = -gstat_data.vol; % reverse contrast for to get sensory drive (passive-fixation)
        end

        % nan out all gstats except those in the frontal lobe label
        inlabel_data = gstat_data.vol(frontal_verts{hh}.Var1+1);
        gstat_data.vol(:) = nan;
        gstat_data.vol(frontal_verts{hh}.Var1+1) = inlabel_data;

        % Normalize gstats
        disp([hemi ' ' contrast ' gstat mean=' num2str(nanmean(gstat_data.vol)) ' max=' num2str(max(gstat_data.vol)) ])
        neg_tstat_mask = gstat_data.vol < 0;
        gstat_data.vol(frontal_verts{hh}.Var1+1) = zscore(gstat_data.vol(frontal_verts{hh}.Var1+1));
        gstat_data.vol(neg_tstat_mask) = nan;
        gstat_maps(:,cc) = gstat_data.vol;

        % nan out all vertices outside union of significant clusters 
        gstat_maps(~union_sig,cc) = nan;

    end

    % Subtract normalized stat maps
    comp_maps(:,hh) = gstat_maps(:,2) - gstat_maps(:,1); 

    % Save as nii
    outfile = ['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/gradient_niftis/' hemi '_' latmed '_' contrasts{2} '-' contrasts{1} '_groupcomp_frontal_tfce.nii.gz'];
    gstat_data.vol = comp_maps(:,hh)';
    MRIwrite(gstat_data, outfile)
end




%% visual - auditory

contrasts = {'aA-aP', 'vA-vP'};
comp_maps = nan(fs_nverts,2);

for hh = 1:2
    hemi = hemis{hh};
    gstat_maps = nan(fs_nverts, 2);

    % Get union of significant tfce clusters for each contrast
    contrast1_tfce_p = MRIread([tfce_path contrasts{1} '/' hemi '/_tfce_tstat_fwep_nifti1.nii']);
    contrast1_tfce_sig = contrast1_tfce_p.vol > 1.3;
    contrast2_tfce_p = MRIread([tfce_path contrasts{2} '/' hemi '/_tfce_tstat_fwep_nifti1.nii']);
    contrast2_tfce_sig = contrast2_tfce_p.vol > 1.3;
    union_sig = contrast1_tfce_sig | contrast2_tfce_sig;

    for cc = 1:2 % contrast
        contrast = contrasts{cc};

        % Load tstats
        gstat_path = [gstat_pathbase hemi '.ces.localizer_groupavg_' contrast '.glmres/osgm/gamma.mgh'];
        gstat_data = MRIread(gstat_path);

        if contains(contrast, 'f-')
            gstat_data.vol = -gstat_data.vol; % reverse contrast for to get sensory drive (passive-fixation)
        end

        % nan out all gstats except those in the frontal lobe label
        inlabel_data = gstat_data.vol(frontal_verts{hh}.Var1+1);
        gstat_data.vol(:) = nan;
        gstat_data.vol(frontal_verts{hh}.Var1+1) = inlabel_data;

        % Normalize gstats
        disp([hemi ' ' contrast ' gstat mean=' num2str(nanmean(gstat_data.vol)) ' max=' num2str(max(gstat_data.vol)) ])
        neg_tstat_mask = gstat_data.vol < 0;
        gstat_data.vol(frontal_verts{hh}.Var1+1) = zscore(gstat_data.vol(frontal_verts{hh}.Var1+1));
        gstat_data.vol(neg_tstat_mask) = nan;
        gstat_maps(:,cc) = gstat_data.vol;

        % nan out all vertices outside union of significant clusters 
        gstat_maps(~union_sig,cc) = nan;

    end

    % Subtract normalized stat maps
    comp_maps(:,hh) = gstat_maps(:,1) - gstat_maps(:,2); 

    % Save as nii
    outfile = ['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/gradient_niftis/' hemi '_' latmed '_' contrasts{1} '-' contrasts{2} '_groupcomp_frontal_tfce.nii.gz'];
    gstat_data.vol = comp_maps(:,hh)';
    MRIwrite(gstat_data, outfile)
end


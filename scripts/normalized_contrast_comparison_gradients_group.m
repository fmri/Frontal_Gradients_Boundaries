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
ROI_path = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/';

hemis = {'lh', 'rh'};
modality = 't'; % visual = v, auditory = a, tactile = t

lobe_annot_lh = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/fsaverage/label/lh.PALS_B12_Lobes.annot';
[~, labels, ctable] = read_annotation(lobe_annot_lh);
frontal_code = ctable.table(ismember(ctable.struct_names, 'LOBE.FRONTAL'), 5);
frontal_verts_lh = find(labels==frontal_code)-1; % annots and labels are off by one so subtract 1
lobe_annot_rh = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/fsaverage/label/rh.PALS_B12_Lobes.annot';
[~, labels, ctable] = read_annotation(lobe_annot_rh);
frontal_code = ctable.table(ismember(ctable.struct_names, 'LOBE.FRONTAL'), 5);
frontal_verts_rh = find(labels==frontal_code)-1; % annots and labels are off by one so subtract 1
frontal_verts = {frontal_verts_lh, frontal_verts_rh};

%% Loop through subjs, get t-stats for each contrast, normalize them, then subtract them and save as nii

contrasts = {['f-' modality 'P'], [modality 'A-' modality 'P']};
comp_maps = nan(fs_nverts,2);

for hh = 1:2
    hemi = hemis{hh};
    gstat_maps = nan(fs_nverts, 2);
    for cc = 1:2 % contrast
        contrast = contrasts{cc};

        % Load tstats
        gstat_path = [gstat_pathbase hemi '.ces.localizer_groupavg_' contrast '.glmres/osgm/gamma.mgh'];
        gstat_data = MRIread(gstat_path);

        if cc == 1
            gstat_data.vol = -gstat_data.vol; % reverse contrast for to get sensory drive (passive-fixation)
        end

        % nan out all gstats except those in the label
        inlabel_data = gstat_data.vol(frontal_verts{hh}+1);
        gstat_data.vol(:) = nan;
        gstat_data.vol(frontal_verts{hh}+1) = inlabel_data;

        % Normalize gstats
        disp([hemi ' ' contrast ' gstat mean=' num2str(nanmean(gstat_data.vol)) ' max=' num2str(max(gstat_data.vol)) ])
        gstat_data.vol(frontal_verts{hh}+1) = zscore(gstat_data.vol(frontal_verts{hh}+1));
        gstat_maps(:,cc) = gstat_data.vol;

        % nan out all vertices outside probabilistic
        if strcmp(modality, 'v')
            frontal_label = readtable([ROI_path hemi '.frontal_vis_WM_sensory_combined.label'], 'FileType','text');
        elseif strcmp(modality, 'a')
            frontal_label = readtable([ROI_path hemi '.frontal_aud_WM_sensory_combined.label'], 'FileType','text');
        elseif strcmp(modality, 't')
            frontal_label = readtable([ROI_path hemi '.frontal_tac_WM_sensory_combined.label'], 'FileType','text');
        else
            error('modality not recognized');
        end

        frontal_label_mask = ismember(1:fs_nverts, frontal_label{:,1}+1);
        gstat_maps(~frontal_label_mask,cc) = nan;

    end

    % Subtract normalized stat maps
    comp_maps(:,hh) = gstat_maps(:,2) - gstat_maps(:,1); 

    % Save as nii
    outfile = ['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/gradient_niftis/' hemi '_' contrasts{2} '-' contrasts{1} '_groupcomp_frontal.nii.gz'];
    gstat_data.vol = comp_maps(:,hh)';
    MRIwrite(gstat_data, outfile)
end




%% Same thing except contrast visual and auditory

contrasts = {'f-vP', 'f-aP'};
comp_maps = nan(fs_nverts,2);

for hh = 1:2
    hemi = hemis{hh};
    gstat_maps = nan(fs_nverts, 2);
    for cc = 1:2 % contrast
        contrast = contrasts{cc};

        % Load tstats
        gstat_path = [gstat_pathbase hemi '.ces.localizer_groupavg_' contrast '.glmres/osgm/gamma.mgh'];
        gstat_data = MRIread(gstat_path);

        if contains(contrast, 'f-')
            gstat_data.vol = -gstat_data.vol; % reverse contrast for to get sensory drive (passive-fixation)
        end

        % nan out all gstats except those in the label
        inlabel_data = gstat_data.vol(frontal_verts{hh}+1);
        gstat_data.vol(:) = nan;
        gstat_data.vol(frontal_verts{hh}+1) = inlabel_data;

        % Normalize gstats
        disp([hemi ' ' contrast ' gstat mean=' num2str(nanmean(gstat_data.vol)) ' max=' num2str(max(gstat_data.vol)) ])
        gstat_data.vol(frontal_verts{hh}+1) = zscore(gstat_data.vol(frontal_verts{hh}+1));
        gstat_maps(:,cc) = gstat_data.vol;

        % nan out all vertices outside probabilistic
        frontal_label = readtable([ROI_path hemi '.frontal_VisAud_sensory_combined.label'], 'FileType','text');
        frontal_label_mask = ismember(1:fs_nverts, frontal_label{:,1}+1);
        gstat_maps(~frontal_label_mask,cc) = nan;

    end

    % Subtract normalized stat maps
    comp_maps(:,hh) = gstat_maps(:,2) - gstat_maps(:,1); 

    % Save as nii
    outfile = ['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/gradient_niftis/' hemi '_' contrasts{2} '-' contrasts{1} '_groupcomp_frontal.nii.gz'];
    gstat_data.vol = comp_maps(:,hh)';
    MRIwrite(gstat_data, outfile)
end



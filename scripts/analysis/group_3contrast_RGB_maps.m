%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to map 3 fMRI contrasts into RGB values
%%% that will display all contrasts' information by vertex/voxel
%%% Tom Possidente - July 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/functions/');
ccc;

%% Initialize Key Variables
contrasts = {'aP-f', 'vAaA-vPaP', 'vP-f'};
region = 'frontal_lateral_cortex'; % "frontal_lateral_cortex", "frontal_medial_cortex", "pVis_probabilistic", "pAud_probabilistic"

hemis = {'lh', 'rh'};
fs_num = 163842;

RGB_data = nan(fs_num, 3, 2); % vertices x contrasts x hemispheres
contrast_stat_data = nan(fs_num,3,2);
contrast_stat_data_cut = nan(fs_num,3,2);

N_contrasts = length(contrasts);

stat_pathbase = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/grouplevel/';
tfce_path = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/group_contrast_TFCE/';

% Load search space region for each hemisphere
region_label_lh = ['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/lh.' region '.label'];
region_verts_lh = readtable(region_label_lh, 'FileType','text');
region_label_rh = ['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/rh.' region '.label'];
region_verts_rh = readtable(region_label_rh, 'FileType','text');
region_verts = {region_verts_lh, region_verts_rh};

%% Load contrast data and convert to RGB range
% Get union of significant tfce clusters for each contrast

for hh = 1:2
    hemi = hemis{hh}; 

    contrast1_tfce_p = MRIread([tfce_path contrasts{1} '/' hemi '/_tfce_tstat_fwep_nifti1.nii']);
    contrast1_tfce_sig = contrast1_tfce_p.vol > 1.3;
    contrast2_tfce_p = MRIread([tfce_path contrasts{2} '/' hemi '/_tfce_tstat_fwep_nifti1.nii']);
    contrast2_tfce_sig = contrast2_tfce_p.vol > 1.3;
    contrast3_tfce_p = MRIread([tfce_path contrasts{3} '/' hemi '/_tfce_tstat_fwep_nifti1.nii']);
    contrast3_tfce_sig = contrast3_tfce_p.vol > 1.3;
    union_sig = contrast1_tfce_sig | contrast2_tfce_sig | contrast3_tfce_sig;

    for cc = 1:3
        contrast = contrasts{cc};
        
        if contains(contrast, '-f') % switch contrast string order bc it's switched in the file name
            contrast_str = split(contrast,'-');
            contrast_str = [contrast_str{2} '-' contrast_str{1}];
            reverse_contrast = true;
        else
            contrast_str = contrast;
            reverse_contrast = false;
        end
        stat_path = [stat_pathbase hemi '.ces.localizer_groupavg_' contrast_str '.glmres/osgm/z.mgh'];
        MRIdata = MRIread(stat_path);
        if reverse_contrast
            MRIdata.vol = -MRIdata.vol;
        end

        % Nan out all vertices not is search space
        inlabel_data = MRIdata.vol(region_verts{hh}.Var1+1);
        MRIdata.vol(:) = nan;
        MRIdata.vol(region_verts{hh}.Var1+1) = inlabel_data;

        contrast_stat_data(:,cc,hh) = MRIdata.vol;

        % Zero out vertices in search space with negative statistic
        neg_stats = MRIdata.vol<0;
        MRIdata.vol(neg_stats) = nan;

        % Zero out vertices in search space that are not significant under any of the 3 contrasts
        MRIdata.vol(~union_sig) = nan;

        contrast_stat_data_cut(:,cc,hh) = MRIdata.vol;

        % Normalize remaining values between 0-255
        ptiles = prctile(MRIdata.vol, [10,20,30,40,50,60,70,80,90]);
        ptiles = [0, ptiles, inf];
        RGB_vals = linspace(0,240,10);
        for vv = 1:length(ptiles)-1
            mask = MRIdata.vol >= ptiles(vv) & MRIdata.vol < ptiles(vv+1);
            RGB_data(mask,cc,hh) = RGB_vals(vv);
        end
    end

    figure;
    subplot(1,3,1);
    histogram(contrast_stat_data(:,1,hh), 'NumBins',50); hold on;
    histogram(contrast_stat_data(:,2,hh), 'NumBins',50);
    histogram(contrast_stat_data(:,3,hh), 'NumBins',50);
    legend(contrasts)
    title([hemi ' non-normalized'])
    subplot(1,3,2);
    histogram(contrast_stat_data_cut(:,1,hh), 'NumBins',50); hold on;
    histogram(contrast_stat_data_cut(:,2,hh), 'NumBins',50);
    histogram(contrast_stat_data_cut(:,3,hh), 'NumBins',50);
    title([hemi ' non-normalized cut'])
    subplot(1,3,3);
    histogram(RGB_data(:,1,hh), 'NumBins',50); hold on;
    histogram(RGB_data(:,2,hh), 'NumBins',50);
    histogram(RGB_data(:,3,hh), 'NumBins',50);
    title([hemi ' normalized'])

end


%% Reorder vertices to match label order
lh_reorder = RGB_data(region_verts_lh.Var1+1,:,1);
rh_reorder = RGB_data(region_verts_rh.Var1+1,:,2);


%% Make annotation file
% Read in reference annot file so we can keep the same structure for our annot files
annotpath = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/fsaverage/label/lh.aparc.annot';
[verts_ref, labels_ref, ctable_ref] = read_annotation(annotpath); % Read in lh.aparc.annot file for reference annotation file structure

rgb_data = {lh_reorder, rh_reorder};
annot_outpath = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/';

for hh=1:2
    hemi = hemis{hh};
    rgb = rgb_data{hh};
    rgb_round = round(rgb, -1); % round each RGB to nearest 10 so there aren't so many unique colors
    unique_colors = unique(rgb_round, "rows");
    N_colors = size(unique_colors,1);

    % Set up new annotation structure
    annot_labels = zeros(length(verts_ref),1); % start with labels as all 0s
    annot_ctable.orig_tab = ctable_ref.orig_tab; % keep same colortable reference
    annot_ctable.struct_names = {ctable_ref.struct_names{1}}; % 'unknown' label is always 1st
    annot_ctable.table = ctable_ref.table(1,:); % copy canonical color RGB and label for 'unknown' label

    % Loop over unique colors and add vertices to annotation file
    for uu = 1:N_colors
        color = unique_colors(uu,:);
        color_label = color(1) + color(2)*(2^8) + color(3)*(2^16); % this is how freesurfer determines color labels for some reason
        annot_ctable.struct_names{uu+1} = num2str(uu); % arbitrary name for this label (doesn't matter)
        annot_ctable.table(uu+1,1:3) = color; % add color to table
        annot_ctable.table(uu+1,5) = color_label; % add color label to same row

        vert_mask = rgb_round(:,1)==color(1) & rgb_round(:,2)==color(2) & rgb_round(:,3)==color(3); % get verts with this color
        vert_inds = region_verts{hh}.Var1(vert_mask) + 1; % get vert inds (add one bc vert inds start at 0 in labels)
        assert(all(annot_labels(vert_inds)==0), 'overlap error');
        annot_labels(vert_inds) = color_label;
    end

    % Set number of entries in annot file
    annot_ctable.numEntries = N_colors + 1; % ROIs plus one for 'unknown' label 0

    % Finally write annotation file
    write_annotation([annot_outpath hemi '.'  'VisAudWM_3contrast_RGB_' region '_sensory.annot'], verts_ref, annot_labels, annot_ctable);
    clear annot_ctable
end







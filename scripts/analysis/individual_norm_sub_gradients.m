%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to take contrasts at the individual level
% and take their normalized difference to see if there is evidence for
% gradients of functional properties
% Tom Possidente - June 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Set key variables
experiment_name = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
%subjCodes = subjDf_cut.subjCode(~ismember(subjDf_cut.subjCode, {'AH', 'SL', 'RR'}));
subjCodes = subjDf_cut.subjCode(~ismember(subjDf_cut.subjCode, {'AH', 'SL', 'RR', 'PP', 'MM'}));
N_subjs = length(subjCodes);

fs_nverts = 163842;

contrast_path = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';

contrasts = {'vAaA-vPaP', 'aPvP-f'};
hemis = {'lh', 'rh'};
comp_maps = nan(fs_nverts, 2, N_subjs);

%% Load frontal cortex labels
% latmed = 'medial';
frontal_label_lh = ['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/lh.frontal_cortex.label'];
frontal_verts_lh = readtable(frontal_label_lh, 'FileType','text');
frontal_label_rh = ['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/rh.frontal_cortex.label'];
frontal_verts_rh = readtable(frontal_label_rh, 'FileType','text');
frontal_verts = {frontal_verts_lh, frontal_verts_rh};


%% Loop through subjs and extract normalized difference of contrasts in frontal cortex
for ss = 1:N_subjs
    subjCode = subjCodes{ss};
    for hh = 1:2
        hemi = hemis{hh};

        contrast_data{1} = MRIread([contrast_path subjCode '/localizer/localizer_contrasts_' hemi '/' contrasts{1} '/sig.nii.gz']);
        if contains(contrasts{1}, 'f-') % switch contrast if it's backward
            contrast_data{1}.vol = -contrast_data{1}.vol;
        end
        contrast1_sig = contrast_data{1}.vol > 1.3;

        contrast_data{2} = MRIread([contrast_path subjCode '/localizer/localizer_contrasts_' hemi '/' contrasts{2} '/sig.nii.gz']);
        if contains(contrasts{2}, 'f-') % switch contrast if it's backward
            contrast_data{2}.vol = -contrast_data{2}.vol;
        end
        contrast2_sig = contrast_data{2}.vol > 1.3;
        
        % Get intersection, contrast1 exclusive, contras2 exclusive
        intersection = contrast1_sig & contrast2_sig;
        contrast1_excl = contrast1_sig & ~contrast2_sig;
        contrast2_excl = contrast2_sig & ~contrast1_sig;
        if sum(intersection) < 100
            disp([subjCode ' ' hemi ' intersection: ' num2str(sum(intersection))])
        end
        if sum(contrast1_excl) < 100
            disp([subjCode ' ' hemi ' vAaA-vPaP: ' num2str(sum(contrast1_excl))])
        end        
        if sum(contrast2_excl) < 100
            disp([subjCode ' ' hemi ' vPaP-f: ' num2str(sum(contrast2_excl))])
        end

        % Create new label files for individual contrast significant regions
        % mask1 = ismember(frontal_verts{hh}{:,1}, find(contrast1_excl)-1); 
        % contrast1_label = frontal_verts{hh}(mask1,:);
        % label_fname = ['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/' hemi '.' subjCode '_' contrasts{1} '_frontal_' latmed  '_sig.label'];
        % label_file = fopen(label_fname,'w');
        % fprintf(label_file, ['#!ascii label  , from subject  vox2ras=TkReg\n' num2str(size(contrast1_label,1)) '\n']);
        % writetable(contrast1_label, label_fname, 'Delimiter', 'tab', 'WriteMode', 'append', 'FileType', 'text');
        % fclose(label_file);
        % 
        % mask2 = ismember(frontal_verts{hh}{:,1}, find(contrast2_excl)-1);
        % contrast2_label = frontal_verts{hh}(mask2,:);
        % label_fname = ['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/' hemi '.' subjCode '_' contrasts{2} '_frontal_' latmed  '_sig.label'];
        % label_file = fopen(label_fname,'w');
        % fprintf(label_file, ['#!ascii label  , from subject  vox2ras=TkReg\n' num2str(size(contrast2_label,1)) '\n']);
        % writetable(contrast2_label, label_fname, 'Delimiter', 'tab', 'WriteMode', 'append', 'FileType', 'text');
        % fclose(label_file);
        % 
        % mask12 = ismember(frontal_verts{hh}{:,1}, find(intersection)-1);
        % contrast12_label = frontal_verts{hh}(mask12,:);
        % label_fname = ['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/' hemi '.' subjCode '_' contrasts{1} '_intersection_' contrasts{2} '_frontal_' latmed  '_sig.label'];
        % label_file = fopen(label_fname,'w');
        % fprintf(label_file, ['#!ascii label  , from subject  vox2ras=TkReg\n' num2str(size(contrast12_label,1)) '\n']);
        % writetable(contrast12_label, label_fname, 'Delimiter', 'tab', 'WriteMode', 'append', 'FileType', 'text');
        % fclose(label_file);

        
        % Get union of significant regions
        union_sig = contrast1_sig | contrast2_sig;

        stat_maps = nan(fs_nverts, 2);

        for cc = 1:2
            contrast = contrasts{cc};

            % Get t-stats
            tstat_path = [contrast_path subjCode '/localizer/localizer_contrasts_' hemi '/' contrasts{cc} '/t.nii.gz'];
            tstat_data = MRIread(tstat_path);

            if contains(contrast, 'f-')
                tstat_data.vol = -tstat_data.vol; % reverse contrast for to get sensory drive (passive-fixation)
            end

            % nan out all gstats except those in the frontal lobe label
            inlabel_data = tstat_data.vol(frontal_verts{hh}.Var1+1);
            tstat_data.vol(:) = nan;
            tstat_data.vol(frontal_verts{hh}.Var1+1) = inlabel_data;

            % Normalize tstats
            disp([hemi ' ' contrast ' tstat mean=' num2str(nanmean(tstat_data.vol)) ' max=' num2str(max(tstat_data.vol)) ])
            neg_tstat_mask = tstat_data.vol < 0;
            tstat_data.vol(frontal_verts{hh}.Var1+1) = zscore(tstat_data.vol(frontal_verts{hh}.Var1+1));
            tstat_data.vol(neg_tstat_mask) = nan;
            stat_maps(:,cc) = tstat_data.vol;

            % nan out all vertices outside union of significant clusters
            stat_maps(~union_sig,cc) = nan;

        end

        % Subtract normalized stat maps
        comp_maps(:,hh,ss) = stat_maps(:,1) - stat_maps(:,2);

        % Save as nii
        outfile = ['/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/gradient_niftis/' subjCode '_' hemi '_' contrasts{1} '-' contrasts{2} '_norm_frontal.nii.gz'];
        contrast_data{cc}.vol = comp_maps(:,hh,ss)';
        MRIwrite(contrast_data{cc}, outfile)

    end

end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to create individual subj labels for the
%%% overlap between significant sensory and WM regions within each ROI.
%%%
%%% Tom Possidente - June 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

%% Initialize Key Variables
hemis = {'lh', 'rh'};
num_hemis = length(hemis);

subjCodes = {'AB', 'AF', 'PQ'};
N_subjs = length(subjCodes);

contrasts = {'aA-aP', 'aP-f'};
keyword = 'auditory';

ROIs = {'preSMA', 'inf_lat_frontal', 'aINS', 'sup_lat_frontal', 'midIFS', 'aIPS'}; 
N_ROIs = length(ROIs);

t_thresh = 2;

ROI_dir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/subj_specific_01/';
data_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/unpacked_data_nii_fs_localizer/';
outdir = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/individual_SMCWM_overlap_perROI/';

dice_coefs = nan(N_subjs, N_ROIs, num_hemis);

%% Loop over subjs/hemis and calculate dice coef between contrasts

for ss = 1:N_subjs
    subjCode = subjCodes{ss};
    for hh = 1:num_hemis
        hemi = hemis{hh};
        for rr = 1:N_ROIs
            ROI_name = ROIs{rr};
            label_file = [ROI_dir hemi '.' ROI_name '_' keyword '_' subjCode '.label'];
            if ~isfile(label_file)
                continue;
            end
            ROI = readtable(label_file, "FileType","text");
            vertices = true(height(ROI),1);
            for cc = 1:2
                if endsWith(contrasts{cc}, '-f')
                    contrast = ['f-' contrasts{cc}(1:2)];
                    tmap = MRIread([data_dir subjCode '/localizer/localizer_contrasts_sm0_' hemi '/' contrast '/t.nii.gz']);
                    tmap.vol = -tmap.vol;
                else
                    tmap = MRIread([data_dir subjCode '/localizer/localizer_contrasts_sm0_' hemi '/' contrasts{cc} '/t.nii.gz']);
                end
                tmap_inROI = tmap.vol(ROI.Var1+1);
                vertices(tmap_inROI<t_thresh) = false;
            end
            disp(sum(vertices))
            label_overlap = table2array(ROI(vertices,:));
            label_fname = [outdir subjCode '_' hemi '_' keyword '_' ROI_name '_WMSMCoverlap.label'];
            label_file = fopen(label_fname,'w');
            fprintf(label_file, ['#!ascii label  , from subject  vox2ras=TkReg\n' num2str(size(label_overlap,1)) '\n']);
            writematrix(label_overlap, label_fname, 'Delimiter', 'tab', 'WriteMode', 'append', 'FileType', 'text');
            fclose(label_file);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to make task contrasts out of the
%%% freesurfer preprocessed localizer data surfaces using half of the data
%%% each time (split-half analysis).
%%%
%%% Tom Possidente - September 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/functions/');
ccc;

%% Set up directories and subj info

projectDir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/';

subjCodes = {'MK','AB', 'AD','LA','AE','TP','NM','AF','AG','AI','GG','UV','PQ','KQ','LN','RT','PT','PL','NS','SL'};
N = length(subjCodes);

fsaverage = true;

keyword = 'even';

data_dir = [projectDir 'data/unpacked_data_nii_fs_localizer/'];
analysis_name_lh = ['localizer_contrasts_0sm_lh_' keyword 'runs'];
analysis_name_rh = ['localizer_contrasts_0sm_rh_' keyword 'runs'];
TR = 2;
design_type = 'blocked';
para_name = 'localizer_condition_timing.para';
if fsaverage
    space = 'fsaverage';
else
    space = 'self';
end
funcstem_lh = ['fmcpr_tu.siemens.sm0.' space '.lh.nii.gz'];
funcstem_rh = ['fmcpr_tu.siemens.sm0.' space '.rh.nii.gz'];
rlf_name = ['localizer_contrasts_runlistfile_' keyword 's.txt'];

run_mkanalysis = true;
run_mkcontrast = true;
use_3way_contrasts = false;

if use_3way_contrasts == true
    nconditions = 5;
else
    nconditions = 7;
end


% localizer condition orders (1way)
% 1 = vA
% 2 = vP
% 3 = aA
% 4 = aP
% 5 = tA
% 6 = tP
% 7 = f

% localizer condition orders (3way)
% 1 = aA
% 2 = tA
% 3 = vA
% 4 = f
% 5 = p


%% Create analysis/contrast info

% run make analysis lh/rh
% using -event-related even
% though this is a block design because this option only matters for
% preprocessing which has already been done. Using polyfit 1 (can use 2 if
% we want more noise regressors). spmhrf 0 to use 0 derivatives of the hrf.
% mcextreg option will use 12 motion regressors (3 translation 3 rotation
% and derivatives)

if run_mkanalysis
    unix(['mkanalysis-sess -a ' analysis_name_lh ' -funcstem ' funcstem_lh ...
        ' -surface ' space ' lh -fsd localizer -event-related -paradigm ' para_name ...
        ' -nconditions ' num2str(nconditions) ' -refeventdur 32 -TR ' num2str(TR) ...
        ' -polyfit 1 -spmhrf 0 -mcextreg -runlistfile ' rlf_name ' -per-run -force'])

    unix(['mkanalysis-sess -a ' analysis_name_rh ' -funcstem ' funcstem_rh ...
        ' -surface ' space ' rh -fsd localizer -event-related -paradigm ' para_name ...
        ' -nconditions ' num2str(nconditions) ' -refeventdur 32 -TR ' num2str(TR) ...
        ' -polyfit 1 -spmhrf 0 -mcextreg -runlistfile ' rlf_name ' -per-run -force'])
end


% run make contrast

% localizer condition orders
% 1 = vA
% 2 = vP
% 3 = aA
% 4 = aP
% 5 = tA
% 6 = tP
% 7 = f


if run_mkcontrast
    if use_3way_contrasts
        error('3way contrasts not set up')
    else

        unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
            'vAvP-f -a 1 -a 2 -c 7'])
        unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
            'vAvP-f -a 1 -a 2 -c 7'])

        unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
            'aAaP-f -a 3 -a 4 -c 7'])
        unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
            'aAaP-f -a 3 -a 4 -c 7'])

        unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
            'aA-aP -a 3 -c 4'])
        unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
            'aA-aP -a 3 -c 4'])

        unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
            'aP-f -a 4 -c 7'])
        unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
            'aP-f -a 4 -c 7'])

        unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
            'vA-vP -a 1 -c 2'])
        unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
            'vA-vP -a 1 -c 2'])

        unix(['mkcontrast-sess -analysis ' analysis_name_lh ' -contrast ' ...
            'vP-f -a 2 -c 7'])
        unix(['mkcontrast-sess -analysis ' analysis_name_rh ' -contrast ' ...
            'vP-f -a 2 -c 7'])

    end

end

%% Loop through subjs
parfor ss = 1:N

    subjCode = subjCodes{ss};

    % run glm
    unix(['selxavg3-sess -s ' subjCode ' -d ' data_dir ' -analysis ' analysis_name_lh...
        ' -no-preproc -overwrite'])
    unix(['selxavg3-sess -s ' subjCode ' -d ' data_dir ' -analysis ' analysis_name_rh...
        ' -no-preproc -overwrite'])

end


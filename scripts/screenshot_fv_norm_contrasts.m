%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to call freeview_screenshots to take
% screenshots of all ROIs on all subjs and save them.
%
% Tom Possidente August 1 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/'))
ccc;

%% Get subjIDs
experiment_name = 'spacetime';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.([experiment_name,'Runs']),''),:);
subjCodes = subjDf_cut.subjCode;
subjCodes = subjDf_cut.subjCode(~ismember(subjDf_cut.subjCode, {'AH', 'SL', 'RR', 'PP', 'MM'}));
N = length(subjCodes);


%% Set ROI and contrast lists

contrast_list = {'vAaA-vPaP-aPvP-f'};


lh_ROI_paths = {};
rh_ROI_paths = {};
contrasts = cell(N,2);

roi_base_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/ROIs/';

all_roi_files = {dir(roi_base_dir).name};

for ss = 1:N % loop through subjs
    subjCode = subjCodes{ss};
    subj_ROIs = all_roi_files(contains(all_roi_files, subjCode) & ~contains(all_roi_files,'indiv')); % find all ROI files for subjs
    lh_ROI_files = subj_ROIs(contains(subj_ROIs, 'lh')); % find all lh ROIs
    rh_ROI_files = subj_ROIs(contains(subj_ROIs, 'rh')); % find all rh ROIs
    
    contrasts{ss,1} = 'A-P';
    if ismember(subjCode, {'RR','MM'})
        contrasts{ss,2} = 'vA-aA';
    else
        contrasts{ss,2} = 'V-A';
    end

    for rr = 1:length(lh_ROI_files) % for each ROI
        ROI_split = split(lh_ROI_files{rr}, '_');
        ROI_name = ROI_split{2}; % isolate name of ROI only 
        if ismember(ROI_name, AP_ROIs) % see if it should be displayed on A-P or V-A contrast
            cc = 1;
        elseif ismember(ROI_name, VA_ROIs)
            cc = 2;
        else
            error('Unrecognized ROI');
        end
        lh_ROI_paths{ss, cc}{rr} = [roi_base_dir lh_ROI_files{rr}]; % put it in correct cell for subj and contrast
    end
    lh_ROI_paths{ss,1} = lh_ROI_paths{ss,1}(~cellfun('isempty',lh_ROI_paths{ss,1})); % delete all empty cells - bc I'm bad at coding and couldn't get it to append nicely without creating empty cells :(
    lh_ROI_paths{ss,2} = lh_ROI_paths{ss,2}(~cellfun('isempty',lh_ROI_paths{ss,2}));

    for rr = 1:length(rh_ROI_files)
        ROI_split = split(rh_ROI_files{rr}, '_');
        ROI_name = ROI_split{2};
        if ismember(ROI_name, AP_ROIs)
            cc = 1;
        elseif ismember(ROI_name, VA_ROIs)
            cc = 2;
        else
            error('Unrecognized ROI');
        end
        rh_ROI_paths{ss, cc}{rr} = [roi_base_dir rh_ROI_files{rr}];
    end
    rh_ROI_paths{ss,1} = rh_ROI_paths{ss,1}(~cellfun('isempty',rh_ROI_paths{ss,1}));
    rh_ROI_paths{ss,2} = rh_ROI_paths{ss,2}(~cellfun('isempty',rh_ROI_paths{ss,2}));

end

%% Actually run freeview_screenshots()

freeview_screenshots(subjCodes, contrasts, lh_ROI_paths, rh_ROI_paths, ctable)
crop_ppt_fv_images(subjCodes, contrasts)

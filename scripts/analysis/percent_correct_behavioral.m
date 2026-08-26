%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The purpose of this script is to load in the behavioral data for the
% localizer task and calculate the % correct for each subj for each
% condition for each modality and run statistical testing on them
%
% Created: Tom Possidente - Jan 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/');
ccc;

save_out = false;

%% Load subject info
taskname = 'x1WayLocalizer';
subjDf = load_subjInfo();
subjDf_cut = subjDf(~strcmp(subjDf.('spacetimeRuns'),''),:);
subjCodes = subjDf_cut.subjCode;
subjCodes = subjCodes(~ismember(subjCodes, {'RR', 'PP', 'MM', 'AH'})); % remove excluded subjs
n = length(subjCodes);

%% Loop through subjs and get % correct
behavioral_dir = '/projectnb/somerslab/tom/projects/sensory_networks_FC/data/behavioral/behavioral/';
modality_names = {'visual', 'auditory', 'tactile'};
n_conditions = 2*length(modality_names);
correct_colnames = {'even_trial_responseCorrect', 'odd_trial_responseCorrect'}; % column names with correct/incorrect info
perc_correct_all = nan(n,n_conditions);
n_cond = nan(n,n_conditions);
perc_correct_byrun = nan(n,2,6);

condition_str = {'active', 'passive'};
condition_order = {'visual_active', 'visual_passive', 'auditory_active', 'auditory_passive', 'tactile_active', 'tactile_passive'};

% loop over subjs
for ss = 1:n

    subjCode = subjCodes{ss};
    row_ind = strcmp(subjCode,subjDf_cut.subjCode);
    date = subjDf_cut.([taskname 'Date']){row_ind};
    runs = subjDf_cut.([taskname, 'Runs']){row_ind};
    if contains(runs, '/') 
        runs = replace(runs, '/', ','); % still take all runs
    end
    runs = str2num(runs);
    num_runs = length(runs);
    subjID = [date upper(subjCode)];

    files = {dir([behavioral_dir subjID]).name};
    if isempty(files)
        disp(['No behavioral data files for subj ' subjCode]);
        continue
    end

    if strcmp(taskname, 'x1WayLocalizer')
        files_csv = files( contains(files, '.csv') & (contains(files, 'trifloc_task_') | contains(files, 'ac_trifloc_task')) ); % Get only csv files in this dir
    else
        error('unrecognized taskname');
    end

    assert(length(files_csv)==num_runs, ['Subj ' subjCode ': Number of runs specified in subjInfo.csv does not match number of runs found in behavioral data files']);

    responses = [];
    modalities = {};
    conditions = {};

    % Loop through each run and collect reponses, modality, and condition data
    len_runs = zeros(num_runs+1,1);
    for rr = 1:num_runs
        behavioral_data = readtable([behavioral_dir subjID, '/', files_csv{rr}], 'Delimiter',',');
        nrows = height(behavioral_data);
        modalities(ss,rr,1,:) = behavioral_data.modality;
        modalities(ss,rr,2,:) = behavioral_data.modality;
        conditions(ss,rr,1,:) = behavioral_data.type;
        conditions(ss,rr,2,:) = behavioral_data.type;
        cue_text_stopped(rr,1,:) = behavioral_data.cue_text_stopped;
        cue_text_stopped(rr,2,:) = behavioral_data.cue_text_stopped;

        which_column = find(ismember(correct_colnames, behavioral_data.Properties.VariableNames));
        responses(ss,rr,1,:) = behavioral_data.(correct_colnames{which_column(1)});
        responses(ss,rr,2,:) = behavioral_data.(correct_colnames{which_column(2)});
    end

    subj_responses = squeeze(responses(ss,:,:,:));

    % Loop through each modality and calculate % correct
    for mm = 1:length(modality_names)
        modality_cond1_mask = squeeze(strcmpi(modality_names{mm}, modalities(ss,:,:,:)) & strcmpi(condition_str{1}, conditions(ss,:,:,:))) & ~strcmpi(cue_text_stopped, 'None');
        n_cond(ss,(mm*2)-1) = sum(modality_cond1_mask,'all');
        perc_correct = mean(subj_responses(modality_cond1_mask), 'all', 'omitnan');
        perc_correct_all(ss,(mm*2)-1) = perc_correct;
        assert(~isnan(perc_correct_all(ss,(mm*2)-1)), 'percent correct calculated as nan, should not happen');

        if perc_correct <= 0.55
            disp(['Subj ' subjCode ' condition ' modality_names{mm} ' cond1 has behavior performance below 55% (' num2str(perc_correct) ')']);
        end

        for rr = 1:num_runs
            perc_correct_byrun(ss,mm,rr) = mean(responses(ss,rr, modality_cond1_mask(rr,:)), 'omitnan');
        end

        modality_cond2_mask = squeeze(strcmpi(modality_names{mm}, modalities(ss,:,:,:)) & strcmpi(condition_str{2}, conditions(ss,:,:,:))) & ~strcmpi(cue_text_stopped, 'None');
        n_cond(ss,(mm*2)) = sum(modality_cond2_mask, 'all');
        perc_correct = mean(subj_responses(modality_cond2_mask), 'all', 'omitnan');
        perc_correct_all(ss,(mm*2)) = perc_correct;
        assert(~isnan(perc_correct_all(ss,(mm*2))) || strcmp(condition_str{2}, 'passive'), 'percent correct calculated as nan, should not happen');
        
        if perc_correct <= 0.55
            disp(['Subj ' subjCode ' condition ' modality_names{mm} ' cond2 has behavior performance below 55% (' num2str(perc_correct) ')']);
        end

    end

end

%% Save out results
perc_correct_all = array2table(perc_correct_all(:,[1,3]), 'VariableNames',condition_order([1,3]), 'RowNames',subjCodes)
if save_out
    save('behavioral_percent_correct_data.mat', 'subjCodes', 'perc_correct_all', 'n_cond');
    writetable(perc_correct_all, 'behavioral_data.csv')
end

missing_all_data = all(ismissing(perc_correct_all),2);
perc_correct_all=perc_correct_all(~missing_all_data,:); % delete any rows with all nans
n = n-sum(missing_all_data);

perc_correct_matrix = table2array(perc_correct_all);


%% Run 2 sample matched t-test
[h,p,ci,t] = ttest(perc_correct_all.visual_active, perc_correct_all.auditory_active)
means = [mean(perc_correct_all.visual_active), mean(perc_correct_all.auditory_active)]

%% Plot swarmplot for visualization
close all;
figure;
swarmchart([repmat(1,n,1); repmat(2,n,1)], [perc_correct_all.visual_active; perc_correct_all.auditory_active], 'XJitterWidth',0.3);
hold on;
boxplot(perc_correct_matrix, 'Labels', {'Visual WM','Auditory WM'});
grid on;
ylim([0,1]);
ylabel('Proportion Correct');




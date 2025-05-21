%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to calculate and plot the RAS coordinates
%%% and overlap of sensory and WM frontal probabilistic ROIs to assess potential for
%%% gradients and overlap
%%%
%%% Tom Possidente - May 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/');
ccc;

%% Initialize Key Variables
hemis = {'lh', 'rh'};

ROIs_pathbase = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/';

contrast_groups = {{'f-vP+f-aP', 'vA-vP+aA-aP'}, {'f-vP+f-aP+f-tP', 'vA-vP+aA-aP+tA-tP'}};
contrast_group_strs = {'auditory and visual', 'auditory, tactile, visual'};
RAS_data = cell(2,length(contrast_groups),3); % hemis x contrast_groups x 3 (sensory, overlap, WM)

%% Loop through hemis/contrasts, get RAS coordinates, and calculate overlap

for hh = 1:length(hemis)
    hemi = hemis{hh};
    for cc = 1:length(contrast_groups)
        contrast_group = contrast_groups{cc};
        labels = cell(length(contrast_groups),1);
        for ii = 1:length(contrast_groups)
            contrast = contrast_group{ii};
            if ii == 1
                path = [ROIs_pathbase hemi '.' contrast '_frontal_probabilistic_thresh4.label'];
            else
                path = [ROIs_pathbase hemi '.' contrast '_frontal_probabilistic_thresh6.label'];
            end
            labels{ii} = readtable(path, 'FileType','text');
        end
        sensory_mask = ~ismember(labels{1}{:,1}, labels{2}{:,1});
        overlap_mask = ismember(labels{1}{:,1}, labels{2}{:,1});
        WM_mask = ~ismember(labels{2}{:,1}, labels{1}{:,1});
        RAS_data{hh,cc,1} = labels{1}{sensory_mask,:};
        RAS_data{hh,cc,2} = labels{1}{overlap_mask,:};
        RAS_data{hh,cc,3} = labels{2}{WM_mask,:};
    end
end

%% Calculate significance using anova
group_data = cell(2,1);
hemis = ["lh", "rh"];
types = ["sensory", "overlap", "WM"];
for aa = 2:3 % only AS axes so we can combine hemispheres
    ax_data = [];
    hemisphere = [];
    type = [];
    for hh = 1:length(hemis)
        hemi = hemis(hh);
        for ii = 1:3 % sensory, overlap, WM
            data_curr = RAS_data{hh,1,ii}(:,1+aa);
            ax_data = [ax_data; data_curr];
            hemisphere = [hemisphere; repmat(hemi, length(data_curr), 1)];
            type = [type; repmat(types(ii), length(data_curr), 1)];
        end
    end
    factors = {hemisphere, type};
    tbl = table(ax_data, hemisphere, type, 'VariableNames', {'axis_data', 'hemisphere', 'type'});
    res = anova(tbl, 'axis_data ~ 1 + hemisphere*type')
    group_data{aa} = groupmeans(res,["hemisphere" "type"]);
    multcompare(res, ["hemisphere" "type"], 'CriticalValueType', 'bonferroni')
    %plotComparisons(res,["type","hemisphere"])
end

%% Swarm plots of sensory, overlap, and WM coords
RAS_str = {'Left-Right Coordinate', 'Posterior-Anterior Coordinate', 'Inferior-Superior Coordinate'};
types_str = {'sensory', 'overlap', 'WM'};

for aa = 2:3
    figure;
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        for cc = 1 %1:length(contrast_groups)
            if hh==1
                xpos = [1,3,5];
            else
                xpos = [2,4,6];
            end
            swarmchart(ones(length(RAS_data{hh,cc,1}(:,1+aa)),1)*xpos(1), RAS_data{hh,cc,1}(:,1+aa), [], [0 0.4470 0.7410], 'filled'); hold on;
            swarmchart(ones(length(RAS_data{hh,cc,2}(:,1+aa)),1)*xpos(2), RAS_data{hh,cc,2}(:,1+aa), [], [0.85 0.325 0.098], 'filled');
            swarmchart(ones(length(RAS_data{hh,cc,3}(:,1+aa)),1)*xpos(3), RAS_data{hh,cc,3}(:,1+aa), [], [0.929 0.694 0.125], 'filled');

            row_mask = ismember(group_data{aa}.hemisphere, hemi) & ismember(group_data{aa}.type, types_str{1});
            e = errorbar(xpos(1), group_data{aa}{row_mask,3}, std(RAS_data{hh,cc,1}(:,1+aa)), 'vertical');
            e.Marker = '.'; e.MarkerSize = 50; e.Color = 'r'; e.CapSize = 20;
            row_mask = ismember(group_data{aa}.hemisphere, hemi) & ismember(group_data{aa}.type, types_str{2});
            e = errorbar(xpos(2), group_data{aa}{row_mask,3}, std(RAS_data{hh,cc,2}(:,1+aa)), 'vertical');
            e.Marker = '.'; e.MarkerSize = 50; e.Color = 'r'; e.CapSize = 20;
            row_mask = ismember(group_data{aa}.hemisphere, hemi) & ismember(group_data{aa}.type, types_str{3});
            e = errorbar(xpos(3), group_data{aa}{row_mask,3}, std(RAS_data{hh,cc,3}(:,1+aa)), 'vertical');
            e.Marker = '.'; e.MarkerSize = 50; e.Color = 'r'; e.CapSize = 20;
            
        end
    end
    xticks([1,2,3,4,5,6]);
    xticklabels({'Sensory lh', 'Sensory rh' 'Overlap lh', 'Overlap rh', 'WM lh', 'WM rh'});
    ylabel(RAS_str{aa});
    set(gca, 'FontSize', 18);

end






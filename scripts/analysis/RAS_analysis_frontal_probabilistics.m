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

%contrast_groups = {{'f-vP+f-aP', 'vA-vP+aA-aP'}, {'f-vP+f-aP+f-tP', 'vA-vP+aA-aP+tA-tP'}};
contrast_groups = {{'f-vP+f-aP', 'vA-vP+aA-aP'}};

contrast_group_strs = {'auditory and visual'};
region = 'precentral_sulcus_gyrus'; % lateral, medial, or PCSG (precentral sulcus/gyrus)
RAS_data = cell(2,length(contrast_groups),3); % hemis x contrast_groups x 3 (sensory, overlap, WM)
RAS_table = table();

%% Loop through hemis/contrasts, get RAS coordinates, and calculate overlap

for hh = 1:length(hemis)
    hemi = hemis{hh};
    for cc = 1:length(contrast_groups)
        contrast_group = contrast_groups{cc};
        labels = cell(length(contrast_groups),1);
        for ii = 1:length(contrast_group)
            contrast = contrast_group{ii};
            if ii == 1
                if strcmp(hemi,'rh') && ismember(region, {'lateral', 'precentral_sulcus_gyrus'})
                    path = [ROIs_pathbase hemi '.' contrast '_frontal_' region '_probabilistic_thresh4_removedmotor.label']; % passive
                else
                    path = [ROIs_pathbase hemi '.' contrast '_frontal_' region '_probabilistic_thresh4.label']; % passive
                end
            else
                path = [ROIs_pathbase hemi '.' contrast '_frontal_' region '_probabilistic_thresh6.label']; % active
            end
            labels{ii} = readtable(path, 'FileType','text');
        end
        sensory_mask = ~ismember(labels{1}{:,1}, labels{2}{:,1});
        overlap_mask = ismember(labels{1}{:,1}, labels{2}{:,1});
        WM_mask = ~ismember(labels{2}{:,1}, labels{1}{:,1});
        RAS_data{hh,cc,1} = labels{1}{sensory_mask,:};
        RAS_data{hh,cc,2} = labels{1}{overlap_mask,:};
        RAS_data{hh,cc,3} = labels{2}{WM_mask,:};
        RAS_table = [ RAS_table; array2table( [labels{1}{sensory_mask,:}, ones(sum(sensory_mask),1), ones(sum(sensory_mask),1)*hh] ) ];
        RAS_table = [ RAS_table; array2table( [labels{1}{overlap_mask,:}, ones(sum(overlap_mask),1)*2, ones(sum(overlap_mask),1)*hh] ) ];
        RAS_table = [ RAS_table; array2table( [labels{2}{WM_mask,:}, ones(sum(WM_mask),1)*3, ones(sum(WM_mask),1)*hh] ) ];
    end
end

RAS_table.Properties.VariableNames = {'vertex', 'R', 'A', 'S', 'Nsubjs', 'type', 'hemisphere'};

%% Calculate significance using anova
group_data = cell(2,1);
hemis = ["lh", "rh"];
types = ["sensory", "overlap", "WM"];
for aa = 1:3 % only AS axes so we can combine hemispheres
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

for aa = 1:3
    figure;
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        for cc = 1 %1:length(contrast_groups)
            if hh==1
                xpos = [1,2,3];
            elseif hh==2
                xpos = [5,6,7];
            else
                xpos = [9,10,11];
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
    xticks([1,2,3,5,6,7,9,10,11]);
    xticklabels({'Sensory lh', 'Overlap lh' 'WM lh', 'Sensory rh', 'Overlap rh', 'WM rh'});
    ylabel(RAS_str{aa});
    set(gca, 'FontSize', 18);

end


%% 2D anterior-posterior inferior-suprerior scatter grouped by sensory, overlap, WM (and hemi)
RAS_data = squeeze(RAS_data);
colors = {'r', 'b', 'g'};
score = {};
for hh = 1:length(hemis)
    hemi = hemis{hh};
    t = [];
    figure;
    for cc = 1:3 % loop over sensory, overlap, WM
        t(end+1) = scatter(RAS_data{hh,cc}(:,3), RAS_data{hh,cc}(:,4), 50, colors{cc});
        hold on;
    end
    set(gca, 'xdir', 'reverse');

    % Plot PCs
    hemi_mask = RAS_table.hemisphere == hh;
    [coef, score{hh}, latent, ~, explained, mu] = pca(RAS_table{hemi_mask,[3,4]});
    scale = 3*sqrt(latent);

    q = quiver(mu(1), mu(2), coef(1,1)*scale(1), coef(2,1)*scale(2), 'Color', 'k', 'LineWidth',2, 'MaxHeadSize', 0.5);
    quiver(mu(1), mu(2), -coef(1,1)*scale(1), -coef(2,1)*scale(2), 'Color', 'k', 'LineWidth',2, 'MaxHeadSize', 0.5)

    quiver(mu(1), mu(2), coef(1,2)*scale(2), coef(2,2)*scale(2), 'Color', 'k', 'LineWidth',2, 'MaxHeadSize', 0.5)
    quiver(mu(1), mu(2), -coef(1,2)*scale(2), -coef(2,2)*scale(2), 'Color', 'k', 'LineWidth',2, 'MaxHeadSize', 0.5)
    
    legend([t(1), t(2), t(3), q], {[hemi ' sensory'], [hemi ' overlap'], [hemi ' WM'], 'PC1, PC2'});
    xlabel('Anterior <-> Posterior');
    ylabel('Inferior <-> Superior');
    set(gca, 'FontSize', 18);
end






%% Run ANOVA with hemisphere and type using PC1 coordinate
group_data = {};
for hh = 1:2
    hemi = hemis{hh};
    hemi_mask = RAS_table.hemisphere == hh;
    
    disp(['PC1, ' hemi ' ANOVA']);
    RAS_table.PC1(hemi_mask) = score{hh}(:,1);
    res = anova(RAS_table, 'PC1 ~ 1 + hemisphere*type')
    group_data{1} = groupmeans(res,["hemisphere" "type"]);
    multcompare(res, ["hemisphere" "type"], 'CriticalValueType', 'bonferroni')

    disp(['PC2, ' hemi ' ANOVA']);
    RAS_table.PC2(hemi_mask) = score{hh}(:,2);
    res = anova(RAS_table, 'PC2 ~ 1 + hemisphere*type')
    group_data{2} = groupmeans(res,["hemisphere" "type"]);
    multcompare(res, ["hemisphere" "type"], 'CriticalValueType', 'bonferroni')
end

%% Plot PC1 Swarmplots
types_str = {'sensory', 'overlap', 'WM'};
colors = [0 0.4470 0.7410; 0.85 0.325 0.098; 0.929 0.694 0.125];

figure;
for hh = 1:length(hemis)
    hemi = hemis{hh};
    for tt = 1: length(types_str)
        if hh==1
            xpos = [1,2,3];
        else
            xpos = [5,6,7];
        end

        mask = RAS_table.hemisphere==hh & RAS_table.type == tt;

        swarmchart(ones(sum(mask),1)*xpos(tt), RAS_table.PC1(mask), [], colors(tt,:), 'filled'); hold on;

        e = errorbar(xpos(tt), group_data{1}.Mean(group_data{1}.type==tt & group_data{1}.hemisphere==hh), std(RAS_table.PC1(mask)), 'vertical');
        e.Marker = '.'; e.MarkerSize = 50; e.Color = 'r'; e.CapSize = 20;

    end
end

xticks([1,2,3,5,6,7]);
xticklabels({'Sensory lh', 'Overlap lh' 'WM lh', 'Sensory rh', 'Overlap rh', 'WM rh'});
ylabel('PC1');
set(gca, 'FontSize', 18);


%% Plot PC2 Swarmplots
types_str = {'sensory', 'overlap', 'WM'};
colors = [0 0.4470 0.7410; 0.85 0.325 0.098; 0.929 0.694 0.125];

figure;
for hh = 1:length(hemis)
    hemi = hemis{hh};
    for tt = 1: length(types_str)
        if hh==1
            xpos = [1,2,3];
        else
            xpos = [5,6,7];
        end

        mask = RAS_table.hemisphere==hh & RAS_table.type == tt;

        swarmchart(ones(sum(mask),1)*xpos(tt), RAS_table.PC2(mask), [], colors(tt,:), 'filled'); hold on;

        e = errorbar(xpos(tt), group_data{2}.Mean(group_data{2}.type==tt & group_data{2}.hemisphere==hh), std(RAS_table.PC2(mask)), 'vertical');
        e.Marker = '.'; e.MarkerSize = 50; e.Color = 'r'; e.CapSize = 20;

    end
end

xticks([1,2,3,5,6,7]);
xticklabels({'Sensory lh', 'Overlap lh' 'WM lh', 'Sensory rh', 'Overlap rh', 'WM rh'});
ylabel('PC2');
set(gca, 'FontSize', 18);


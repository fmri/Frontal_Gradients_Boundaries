%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to calculate and plot the RAS coordinates
%%% and overlap of sensory and WM frontal probabilistic ROIs to assess potential for
%%% gradients and overlap
%%%
%%% Tom Possidente - August 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/projectnb/somerslab/tom/projects/sensory_networks_FC/functions/');
ccc;

%% Initialize Key Variables
hemis = {'lh', 'rh'};

ROIs_pathbase = '/projectnb/somerslab/tom/projects/Frontal_Gradients_Boundaries/data/ROIs/';

%contrast_groups = {{'f-vP+f-aP', 'vA-vP+aA-aP'}, {'f-vP+f-aP+f-tP', 'vA-vP+aA-aP+tA-tP'}};
contrast_groups = {{'f-vP+f-aP', 'vA-vP+aA-aP'}};

contrast_group_strs = {'auditory and visual'};
region = 'medial'; % lateral, medial, or PCSG (precentral sulcus/gyrus)
XY_table = table();

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
        sensory_inds = labels{1}{sensory_mask,1};
        overlap_mask = ismember(labels{1}{:,1}, labels{2}{:,1});
        overlap_inds = labels{1}{overlap_mask,1};
        WM_mask = ~ismember(labels{2}{:,1}, labels{1}{:,1});
        WM_inds = labels{2}{WM_mask,1};
        
        flat_patch = read_patch(['/projectnb/somerslab/tom/projects/sensory_networks_FC/data/recons/fsaverage/surf/' hemi '.' region '_VisAudWM_combined_TFCE_flat.patch']);
        sensory_flatpatch_mask = ismember(flat_patch.ind, sensory_inds);
        overlap_flatpatch_mask = ismember(flat_patch.ind, overlap_inds);
        WM_flatpatch_mask = ismember(flat_patch.ind, WM_inds);

        XY_table = [XY_table; array2table([flat_patch.ind(sensory_flatpatch_mask); flat_patch.x(sensory_flatpatch_mask); flat_patch.y(sensory_flatpatch_mask);...
                                          ones(1,sum(sensory_flatpatch_mask)); ones(1,sum(sensory_flatpatch_mask))*hh]' ) ];
        XY_table = [XY_table; array2table([flat_patch.ind(overlap_flatpatch_mask); flat_patch.x(overlap_flatpatch_mask); flat_patch.y(overlap_flatpatch_mask);...
                                          ones(1,sum(overlap_flatpatch_mask))*2; ones(1,sum(overlap_flatpatch_mask))*hh]' ) ];
        XY_table = [XY_table; array2table([flat_patch.ind(WM_flatpatch_mask); flat_patch.x(WM_flatpatch_mask); flat_patch.y(WM_flatpatch_mask);...
                                          ones(1,sum(WM_flatpatch_mask))*3; ones(1,sum(WM_flatpatch_mask))*hh]' ) ];
    end
end

%RAS_table.Properties.VariableNames = {'vertex', 'R', 'A', 'S', 'Nsubjs', 'type', 'hemisphere'};
XY_table.Properties.VariableNames = {'vertex', 'X', 'Y', 'type', 'hemisphere'};

%% Scatter plot
for hh = 1:length(hemis)
    figure;
    tbl = XY_table(XY_table.hemisphere==hh,:);
    scatter(tbl.X(tbl.type==1), tbl.Y(tbl.type==1), 'r')
    hold on;
    scatter(tbl.X(tbl.type==2), tbl.Y(tbl.type==2), 'b')
    scatter(tbl.X(tbl.type==3), tbl.Y(tbl.type==3), 'g')
    legend({'sensory' 'overlap' 'WM'});
    title([hemis{hh} ' ' region]);
    xlabel('Flatmap X Coordinate');
    ylabel('Flatmap Y Coordinate');
end


% 
% % Calculate significance using anova
% group_data = cell(2,1);
% hemis = ["lh", "rh"];
% types = ["sensory", "overlap", "WM"];
% for aa = 1:3 % only AS axes so we can combine hemispheres
%     ax_data = [];
%     hemisphere = [];
%     type = [];
%     for hh = 1:length(hemis)
%         hemi = hemis(hh);
%         for ii = 1:3 % sensory, overlap, WM
%             data_curr = XY_data{hh,1,ii}(:,1+aa);
%             ax_data = [ax_data; data_curr];
%             hemisphere = [hemisphere; repmat(hemi, length(data_curr), 1)];
%             type = [type; repmat(types(ii), length(data_curr), 1)];
%         end
%     end
%     factors = {hemisphere, type};
%     tbl = table(ax_data, hemisphere, type, 'VariableNames', {'axis_data', 'hemisphere', 'type'});
%     res = anova(tbl, 'axis_data ~ 1 + hemisphere*type')
%     group_data{aa} = groupmeans(res,["hemisphere" "type"]);
%     multcompare(res, ["hemisphere" "type"], 'CriticalValueType', 'bonferroni')
%     plotComparisons(res,["type","hemisphere"])
% end
% 
% % Swarm plots of sensory, overlap, and WM coords
% XY_str = {'X', 'Y'};
% types_str = {'sensory', 'overlap', 'WM'};
% 
% for aa = 1:3
%     figure;
%     for hh = 1:length(hemis)
%         hemi = hemis{hh};
%         for cc = 1 %1:length(contrast_groups)
%             if hh==1
%                 xpos = [1,2,3];
%             elseif hh==2
%                 xpos = [5,6,7];
%             else
%                 xpos = [9,10,11];
%             end
% 
%             swarmchart(ones(length(XY_data{hh,cc,1}(:,1+aa)),1)*xpos(1), XY_data{hh,cc,1}(:,1+aa), [], [0 0.4470 0.7410], 'filled'); hold on;
%             swarmchart(ones(length(XY_data{hh,cc,2}(:,1+aa)),1)*xpos(2), XY_data{hh,cc,2}(:,1+aa), [], [0.85 0.325 0.098], 'filled');
%             swarmchart(ones(length(XY_data{hh,cc,3}(:,1+aa)),1)*xpos(3), XY_data{hh,cc,3}(:,1+aa), [], [0.929 0.694 0.125], 'filled');
% 
%             row_mask = ismember(group_data{aa}.hemisphere, hemi) & ismember(group_data{aa}.type, types_str{1});
%             e = errorbar(xpos(1), group_data{aa}{row_mask,3}, std(XY_data{hh,cc,1}(:,1+aa)), 'vertical');
%             e.Marker = '.'; e.MarkerSize = 50; e.Color = 'r'; e.CapSize = 20;
%             row_mask = ismember(group_data{aa}.hemisphere, hemi) & ismember(group_data{aa}.type, types_str{2});
%             e = errorbar(xpos(2), group_data{aa}{row_mask,3}, std(XY_data{hh,cc,2}(:,1+aa)), 'vertical');
%             e.Marker = '.'; e.MarkerSize = 50; e.Color = 'r'; e.CapSize = 20;
%             row_mask = ismember(group_data{aa}.hemisphere, hemi) & ismember(group_data{aa}.type, types_str{3});
%             e = errorbar(xpos(3), group_data{aa}{row_mask,3}, std(XY_data{hh,cc,3}(:,1+aa)), 'vertical');
%             e.Marker = '.'; e.MarkerSize = 50; e.Color = 'r'; e.CapSize = 20;
% 
%         end
%     end
%     xticks([1,2,3,5,6,7,9,10,11]);
%     xticklabels({'Sensory lh', 'Overlap lh' 'WM lh', 'Sensory rh', 'Overlap rh', 'WM rh'});
%     ylabel(RAS_str{aa});
%     set(gca, 'FontSize', 18);
% 
% end
% 
% 
% % 2D anterior-posterior inferior-suprerior scatter grouped by sensory, overlap, WM (and hemi)
% XY_data = squeeze(XY_data);
% colors = {'r', 'b', 'g'};
% score = {};
% for hh = 1:length(hemis)
%     hemi = hemis{hh};
%     t = [];
%     figure;
%     for cc = 1:3 % loop over sensory, overlap, WM
%         t(end+1) = scatter(XY_data{hh,cc}(:,3), XY_data{hh,cc}(:,4), 50, colors{cc});
%         hold on;
%     end
%     set(gca, 'xdir', 'reverse');
% 
%     Plot PCs
%     hemi_mask = XY_table.hemisphere == hh;
%     [coef, score{hh}, latent, ~, explained, mu] = pca(XY_table{hemi_mask,[3,4]});
%     scale = 3*sqrt(latent);
% 
%     q = quiver(mu(1), mu(2), coef(1,1)*scale(1), coef(2,1)*scale(2), 'Color', 'k', 'LineWidth',2, 'MaxHeadSize', 0.5);
%     quiver(mu(1), mu(2), -coef(1,1)*scale(1), -coef(2,1)*scale(2), 'Color', 'k', 'LineWidth',2, 'MaxHeadSize', 0.5)
% 
%     quiver(mu(1), mu(2), coef(1,2)*scale(2), coef(2,2)*scale(2), 'Color', 'k', 'LineWidth',2, 'MaxHeadSize', 0.5)
%     quiver(mu(1), mu(2), -coef(1,2)*scale(2), -coef(2,2)*scale(2), 'Color', 'k', 'LineWidth',2, 'MaxHeadSize', 0.5)
% 
%     legend([t(1), t(2), t(3), q], {[hemi ' sensory'], [hemi ' overlap'], [hemi ' WM'], 'PC1, PC2'});
%     xlabel('Anterior <-> Posterior');
%     ylabel('Inferior <-> Superior');
%     set(gca, 'FontSize', 18);
% end
% 
% 
% 
% 
% 
% 
% % Run ANOVA with hemisphere and type using PC1 coordinate
% group_data = {};
% for hh = 1:2
%     hemi = hemis{hh};
%     hemi_mask = XY_table.hemisphere == hh;
% 
%     disp(['PC1, ' hemi ' ANOVA']);'Inferior-Superior Coordinate'
%     XY_table.PC1(hemi_mask) = score{hh}(:,1);
%     res = anova(XY_table, 'PC1 ~ 1 + hemisphere*type')
%     group_data{1} = groupmeans(res,["hemisphere" "type"]);
%     multcompare(res, ["hemisphere" "type"], 'CriticalValueType', 'bonferroni')
% 
%     disp(['PC2, ' hemi ' ANOVA']);
%     XY_table.PC2(hemi_mask) = score{hh}(:,2);
%     res = anova(XY_table, 'PC2 ~ 1 + hemisphere*type')
%     group_data{2} = groupmeans(res,["hemisphere" "type"]);
%     multcompare(res, ["hemisphere" "type"], 'CriticalValueType', 'bonferroni')
% end
% 
% % Plot PC1 Swarmplots
% types_str = {'sensory', 'overlap', 'WM'};
% colors = [0 0.4470 0.7410; 0.85 0.325 0.098; 0.929 0.694 0.125];
% 
% figure;
% for hh = 1:length(hemis)
%     hemi = hemis{hh};
%     for tt = 1: length(types_str)
%         if hh==1
%             xpos = [1,2,3];
%         else
%             xpos = [5,6,7];
%         end
% 
%         mask = XY_table.hemisphere==hh & XY_table.type == tt;
% 
%         swarmchart(ones(sum(mask),1)*xpos(tt), XY_table.PC1(mask), [], colors(tt,:), 'filled'); hold on;
% 
%         e = errorbar(xpos(tt), group_data{1}.Mean(group_data{1}.type==tt & group_data{1}.hemisphere==hh), std(XY_table.PC1(mask)), 'vertical');
%         e.Marker = '.'; e.MarkerSize = 50; e.Color = 'r'; e.CapSize = 20;
% 
%     end
% end
% 
% xticks([1,2,3,5,6,7]);
% xticklabels({'Sensory lh', 'Overlap lh' 'WM lh', 'Sensory rh', 'Overlap rh', 'WM rh'});
% ylabel('PC1');
% set(gca, 'FontSize', 18);
% 
% 
% % Plot PC2 Swarmplots
% types_str = {'sensory', 'overlap', 'WM'};
% colors = [0 0.4470 0.7410; 0.85 0.325 0.098; 0.929 0.694 0.125];
% 
% figure;
% for hh = 1:length(hemis)
%     hemi = hemis{hh};
%     for tt = 1: length(types_str)
%         if hh==1
%             xpos = [1,2,3];
%         else
%             xpos = [5,6,7];
%         end
% 
%         mask = XY_table.hemisphere==hh & XY_table.type == tt;
% 
%         swarmchart(ones(sum(mask),1)*xpos(tt), XY_table.PC2(mask), [], colors(tt,:), 'filled'); hold on;
% 
%         e = errorbar(xpos(tt), group_data{2}.Mean(group_data{2}.type==tt & group_data{2}.hemisphere==hh), std(XY_table.PC2(mask)), 'vertical');
%         e.Marker = '.'; e.MarkerSize = 50; e.Color = 'r'; e.CapSize = 20;
% 
%     end
% end
% 
% xticks([1,2,3,5,6,7]);
% xticklabels({'Sensory lh', 'Overlap lh' 'WM lh', 'Sensory rh', 'Overlap rh', 'WM rh'});
% ylabel('PC2');
% set(gca, 'FontSize', 18);
% 

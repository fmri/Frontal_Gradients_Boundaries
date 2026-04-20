function [rotated_xs,rotated_ys] = group_find_axis(method, group_data_diff, xs, ys, hemi)
%GROUP_FIND_AXIS attempts to calculate the axis of greatest change in a 2D flat patch 

    arguments (Input)
        method {mustBeText, mustBeNonempty} % which method to use: regression (uses regression to find axis of largest difference) or average (uses weighted average of positive and negative pts to make line)
        group_data_diff (:,1) {mustBeMatrix, mustBeNonNan, mustBeNonempty} % data values corresponding to x and y coordinates
        xs (:,1) {mustBeMatrix, mustBeNonNan, mustBeNonempty} % x coordinates
        ys (:,1) {mustBeMatrix, mustBeNonNan, mustBeNonempty} % y coordinates
        hemi {mustBeText} = '' % which hemisphere is this? Optional, for labeling plots
    end
    
    assert(isequal(size(group_data_diff), size(xs), size(ys)), 'group_data_diff, xs, and ys must be the same size');

    % Precompute 100x100 mesh grid for plotting
    [xq,yq] = meshgrid(linspace(min(xs), max(xs), 100), linspace(min(ys), max(ys), 100));

    % Compute axis of greatest change with selected method
    switch method 
        case 'regression'
            % Linear regression to find 2D plane that best fits 3D data
            X = [xs, ys, ones(length(xs),1)];
            coefs = X \ group_data_diff;
            coefs_norm = coefs(1:2) / norm(coefs(1:2)); % normalize coefficients

            % Plot 3D data surface
            zq = griddata(xs, ys, group_data_diff, xq, yq, 'linear'); % interpolate to 100x100 grid
            figure;
            surf(xq, yq, zq); % plot surface
            hold on;

            % Plot 2D plane of best fit
            zgrid = coefs(1) * xq + coefs(2) * yq + coefs(3); % Get z coord for each xy pair
            mesh(xq,yq,zgrid,'FaceColor', 'g', 'FaceAlpha',0.5, 'EdgeColor','none'); % plot plane

            % Plot 2D line along which there is the greatest change in z
            x0 = mean(xs); % just pick mean x and y coords as middle of line
            y0 = mean(ys);
            z0 = coefs(1)*x0 + coefs(2)*y0 + coefs(3);
            x_pts = linspace(min(xs)-x0, max(xs)-x0, 100); 
            y_pts = linspace(min(ys)-y0, max(ys)-y0, 100);
            x_line = x0 + x_pts * coefs_norm(1);
            y_line = y0 + y_pts * coefs_norm(2);
            plot(x_line, y_line, 'r--', 'LineWidth',2);
            xlabel('x');
            ylabel('y');
            zlabel('T-stat Difference');
            title([hemi ' flatmap with axis plotted']);

            % Get new x and y axes directions where x axis is now 2d line of greatest change in z
            x_ax_new = coefs_norm;
            y_ax_new = [-coefs_norm(2); coefs_norm(1)]; % perpendicular to new_x (rotate 90 deg)

        case 'average'
            % Get weighted average of x,y coords with positive z
            pos_mask = group_data_diff >= 0;
            weighted_pos(1) = sum(xs(pos_mask) .* group_data_diff(pos_mask)) / sum(group_data_diff(pos_mask));
            weighted_pos(2) = sum(ys(pos_mask) .* group_data_diff(pos_mask)) / sum(group_data_diff(pos_mask));

            % Get weighted average of x,y coords with negative z
            neg_mask = group_data_diff < 0;
            weighted_neg(1) = sum(xs(neg_mask) .* group_data_diff(neg_mask)) / sum(group_data_diff(neg_mask));
            weighted_neg(2) = sum(ys(neg_mask) .* group_data_diff(neg_mask)) / sum(group_data_diff(neg_mask));

            % Direction vector from pos to neg means will be new x axis
            v = weighted_pos - weighted_neg;
            x_ax_new = v' / norm(v);  % normalize

            % Perpendicular vector (90 deg rotation)
            y_ax_new = [-x_ax_new(2); x_ax_new(1)];

            % Plot 3D data surface
            figure;
            scatter3(xs, ys, group_data_diff);
            hold on;
            plot([weighted_pos(1), weighted_neg(1)], [weighted_pos(2), weighted_neg(2)], 'LineWidth',5, 'Color', 'r');
            xlabel('x'); ylabel('y'); zlabel('T-stats');
            title([hemi 'flatmap with axis plotted']);

        otherwise
            error('method not recognized')
    end

    % Rotate coordinates so that x axis is axis of greatest change
    xy = [xs'; ys'];
    R = [x_ax_new, y_ax_new]'; % rotation matrix
    xy_rotated = R * xy;
    rotated_xs = xy_rotated(1,:); % new coords
    rotated_ys = xy_rotated(2,:); % new coords

    figure;
    scatter(rotated_xs, rotated_ys, [], group_data_diff);
    clim_set = max(abs(prctile(group_data_diff, [10,90])));
    xlabel('x'); ylabel('y'); cb = colorbar; colormap(redbluedark); clim([-clim_set, clim_set]); ylabel(cb, 'T-stat', 'rotation', 270);
    title([hemi ' flatmap rotated x greatest change'])
end
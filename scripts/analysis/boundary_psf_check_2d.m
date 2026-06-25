%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to convolve an underlying functional
%%% boundary with an estimated fMRI PSF, downsample to voxel resolution,
%%% and fit a hinge model. This will help determine the threshold for
%%% gradient distance required to be confident in a gradient vs boundary. 
%%%
%%% Arthur Sangil Lee - June 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

FWHM = 3.5;
num_runs = 1000; % number of simulations
sigma = FWHM/(2*sqrt(2*log(2)));

% Initialize output array
ramp_lengths = zeros(num_runs, 1);

% Step 1: Define the grid parameters (using vectors to save memory)
finegrid_res = 0.01;
x_vec = single(0:finegrid_res:30); % 30mm by 30mm space
y_vec = single(0:finegrid_res:30)'; % Column vector for implicit expansion

% Downsampling parameters for Step 5
block_size = round(2.2 / finegrid_res); % 2200 pixels
num_blocks_x = floor(length(x_vec) / block_size);
num_blocks_y = floor(length(y_vec) / block_size);

% Anonymous function for the Hinge model (Step 7)
% params = [Lower Bound, Upper Bound, Midpoint, Ramp Length]
hinge_model = @(p, x) p(1) + (p(2) - p(1)) .* max(0, min(1, (x - (p(3) - abs(p(4))/2)) ./ abs(p(4))));

for run = 1:num_runs
    fprintf('Starting simulation run %d of %d...\n', run, num_runs);

    % Step 2: Construct two orthogonal vectors crossing between 10 and 20
    x_c = 10 + 10 * rand();
    y_c = 10 + 10 * rand();

    % Random rotation angle for the vectors
    theta = 2 * pi * rand();

    % Vector 1 (Boundary direction) and Vector 2 (Normal direction)
    v_normal = [-sin(theta), cos(theta)];

    % Step 3: Pick one vector as boundary.
    % We determine sides using the dot product with the normal vector.
    % Using implicit expansion to calculate distance across the whole grid efficiently.
    dist_from_boundary = (x_vec - x_c).*v_normal(1) + (y_vec - y_c).*v_normal(2);

    % Assign weights of 0 and 1
    grid_data = single(dist_from_boundary > 0);

    % if you want to plot it: 
    % subplot(2,2,1)
    % surf(grid_data,'LineStyle','none')

    % Step 4: Smooth with 2D Gaussian
    % Note: sigma for imgaussfilt is in pixels.
    sigma_pixels = sigma / finegrid_res;
    smoothed_data = imgaussfilt(grid_data, sigma_pixels);
    % if you want to plot it: 
    % subplot(2,2,2)
    % surf(smoothed_data,'LineStyle','none')

    % Step 5: Downsample into 2.2 x 2.2 grids via block averaging
    % We crop the excess edges that don't fit perfectly into the 2.2 chunks
    crop_data = smoothed_data(1:num_blocks_y*block_size, 1:num_blocks_x*block_size);

    % Fast block averaging using reshape and mean
    downsampled_data = squeeze(mean(mean(reshape(crop_data, block_size, num_blocks_y, block_size, num_blocks_x), 1), 3));
    % subplot(2,2,3)
    % surf(downsampled_data,'LineStyle','none')

    % Coordinates for the centers of the downsampled coarse grid
    x_coarse = (0.5:num_blocks_x) * 2.2;
    y_coarse = (0.5:num_blocks_y) * 2.2;
    [X_coarse, Y_coarse] = meshgrid(x_coarse, y_coarse);

    % Step 6: Collapse the data into 1D along the normal vector
    % Calculate distance of each coarse pixel to the crossing point along the normal
    dist_1D = (X_coarse(:) - x_c) * v_normal(1) + (Y_coarse(:) - y_c) * v_normal(2);
    vals_1D = downsampled_data(:);

    % Ensure 0s are on the left (negative dist) and 1s are on the right (positive dist)
    if mean(vals_1D(dist_1D > 0)) < mean(vals_1D(dist_1D < 0))
        dist_1D = -dist_1D;
    end

    % Step 7: Fit a hinge function to the 1D data
    % Objective function: Sum of Squared Errors (SSE)
    sse_obj = @(params) sum((vals_1D - hinge_model(params, dist_1D)).^2);

    % Initial guess: [Lower=0, Upper=1, Midpoint=0, Ramp Length=Estimated]
    initial_guess = [0, 1, 0, 3 * sigma];

    % Optimize parameters using fminsearch
    options = optimset('Display', 'off');
    best_params = fminsearch(sse_obj, initial_guess, options);

    % Step 8: Extract the x-length for the linear ramp
    ramp_length = abs(best_params(4));
    ramp_lengths(run) = ramp_length;

    % if you want to plot this
    % subplot(2,2,4)
    % plot(dist_1D, vals_1D, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 4, 'DisplayName', '2.2 Voxel Data'); % Plot 2.2mm voxel data points (blue circles)
    % hold on;
    % x_fit = linspace(min(dist_1D), max(dist_1D), 1000); % Generate smooth x values for the fitted hinge line
    % y_fit = hinge_model(best_params, x_fit); % Generate smooth y values for the fitted hinge line
    % plot(x_fit, y_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Fitted Hinge'); % Plot fitted hinge function (red line)
    % xline(0, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Boundary (Crossing Point)'); % Plot boundary (black line) at distance = 0

    fprintf('Run %d completed. Extracted Ramp Length: %.4f\n', run, ramp_length);
end

fprintf('\nAll simulations completed.\n');

P_95 = prctile(ramp_lengths,95)
P_99 = prctile(ramp_lengths,99)
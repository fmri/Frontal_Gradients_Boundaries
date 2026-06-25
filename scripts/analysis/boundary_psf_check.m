%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The purpose of this script is to convolve an underlying functional
%%% boundary with an estimated fMRI PSF, downsample to voxel resolution,
%%% and fit a hinge model. This will help determine the threshold for
%%% gradient distance required to be confident in a gradient vs boundary. 
%%%
%%% Tom Possidente - June 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/projectnb/somerslab/tom/functions/'));
ccc;

% Set key parameters
plot_diagnostics = true; 
voxel_size = 2.2; %mm
FWHM = 3.5; %mm 3.9 from 3T paper, 3.5 from 1.5T paper, 2mm from 7T paper
sigma = FWHM/(2*sqrt(2*log(2)));
underlying_data_step = 0.001; % mm
x = 0:underlying_data_step:20-underlying_data_step; %mm
T = [zeros(length(x)/2,1); ones(length(x)/2,1)]; 

% Build PSF
kernel_range = -ceil(4*sigma/underlying_data_step) : ceil(4*sigma/underlying_data_step);  % 4 sigma kernel 
kernel_x = kernel_range * underlying_data_step;
psf = exp(-(kernel_x.^2) / (2*sigma^2));
psf = psf / sum(psf);  % normalize so it sums to 1

% convolve boundary with PSF
convolved = conv(T, psf, 'same');

% Loop vairables/storage
sample_offsets = 0:10:2200; % offset starting point of downsampling by 0.1mm each time, up to 2.2mm (voxel)
gradient_lengths = nan(length(sample_offsets),1);

for ss = 1:length(sample_offsets)
    % Sample underlying data with voxel resolution 
    samples = 1+sample_offsets(ss):voxel_size/underlying_data_step:length(x);
    sampled_conv = convolved(samples);
    sampled_x = x(samples);
    
    % Run the hinge model fit on downsampled data
    model = fittype('a*(x < x1) + b*(x > (x1+x2)) + ( a + ( (b-a)/((x1+x2)-x1)  * (x-x1) ) ) * ( (x > x1) & (x < (x1+x2)) )', ... % piecewise function, constant from -Inf to x1, linear from x1 to x2, constant from x2 to Inf. Must include y in the function even if it has no mathematical effect
                'dependent', 'z',...
                'independent', {'x'}, ... % y is not actually used in the function but include it so that we can use all 3D data to fit (not just the x and z coordinates)
                'coefficients', {'a','b','x1', 'x2'}); % a is z from -Inf to x1, b is z from x2 to Inf, x1 is where to start linear piece, x2 is where to end linear piece
    
    startpoint_guesses = [-1.5, 1.5, prctile(x, 25), prctile(x, 50)]; % if no group data, hardcoded guesses
    
    lower_bounds = [-10, -10, min(x), 0];
    upper_bounds = [10, 10, max(x), max(x)-min(x)];
    
    [fit_res_curr, gof_curr, info_curr] = fit(sampled_x', sampled_conv, model, 'StartPoint', startpoint_guesses,...
        'Lower', lower_bounds,...
        'Upper', upper_bounds); 
    parameters = [fit_res_curr.x1, fit_res_curr.x2, fit_res_curr.a, fit_res_curr.b]; 
    
    % Get gradient length
    gradient_lengths(ss) = (parameters(2)+parameters(1)) - parameters(1);

    if plot_diagnostics
        if length(findobj('type','figure')) > 10
            disp('too many figures, not plotting more');
        else
            figure;
            p1 = plot(x, T); 
            hold on;
            p2 = plot(x, convolved);
            p3 = scatter(sampled_x, sampled_conv, 100, 'filled');
            %p4 = scatter(upsampled_x, upsampled_conv, 10, 'filled')
            p4 = plot(fit_res_curr);
            legend({'Underlying Neural Data Boundary', 'Convoled with PSF', 'Sampled Points', 'Hinge Model Fit'});
            title(['Gradient Length: ' num2str(gradient_lengths(ss)) 'mm']);
            grid on;
        end
    end
end

disp(['Mean gradient length: ' num2str(mean(gradient_lengths)) 'mm']);
disp(['95th ptile gradient length = ' num2str(prctile(gradient_lengths, 95)) 'mm']);
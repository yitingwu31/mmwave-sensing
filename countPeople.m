% NChirps => Doppler
% NSamples => Range
% NReceivers => Angle
% Wrapper Function

% ==============================================
% Main Function
% ==============================================
% ------------------------------------------------------------------
% Load DataCube and Parameters
% ------------------------------------------------------------------
% Description: Load dataCube file and output number of people(cluster)
%
% Input:  DataCube file name
% Display: Number of cluster
% ------------------------------------------------------------------
function countPeople (filename)

    % Import radarcube
    data = load(filename); % cells are frames, each cell is (NChirps, NReceivers, NSamples)
    
    % Bin ranges
    radarCube = data.radarCube;
    
    numRangeBins = radarCube.rfParams.numRangeBins;
    numChirps = radarCube.dim.numChirps;
    rangeResolutionsInMeters = radarCube.rfParams.rangeResolutionsInMeters;
    dopplerResolutionMps = radarCube.rfParams.dopplerResolutionMps;
    
    ran_bin = (1:numRangeBins) * rangeResolutionsInMeters;
    % dop_bin = (1:d_bin) * dopplerResolutionMps;
    tot_dop_bin = ((numChirps/2) - (1:numChirps)) * dopplerResolutionMps;
    
    % Adjust max range and doppler
    max_ran = 5;
    max_dop = 3;
    range_indices = ran_bin < max_ran;
    ran_bin_mod = ran_bin(range_indices);
    doppler_indices = (tot_dop_bin >= -1*max_dop) & (tot_dop_bin <= max_dop);
    tot_dop_bin_mod = tot_dop_bin(doppler_indices);
    dop_bin_first = find(doppler_indices, 1, 'first');
    dop_bin_last = find(doppler_indices, 1, 'last');
    NRange = length(ran_bin_mod);
    NDoppler = length(tot_dop_bin_mod);
    
    % Apply windowing and FFT processing
    Nframes = length(radarCube.data);
    [NChirps, NReceivers, NSamples] = size(radarCube.data{1});
    fprintf("NChirps: %d, NReceivers: %d, NSamples: %d\n", NChirps, NReceivers, NSamples);
    
    %Process each frame
    winrange = hann(NSamples);
    winrange3D = permute(repmat(winrange, [1, NChirps, NReceivers]), [2, 1, 3]);
    windoppler = hann(NChirps);
    windoppler3D = repmat(windoppler, [1, NSamples, NReceivers]);
    
    % Individual frames to save
    frame1 = zeros(NChirps, NSamples);
    prevframe = zeros(NChirps, NSamples);
    midframe = zeros(NChirps, NSamples);
    
    % All frames to save
    frames = zeros(Nframes, NChirps, NSamples);
    frames_nozero = zeros(Nframes, NChirps, NSamples);
    
    % figure(1);
    for i=1:Nframes
        frame = permute(radarCube.data{i}, [1, 3, 2]); % permute to match HW3 order of dimensions which is (NChirps, NSamples, NReceivers)
        % undo range FFT
        frame = ifft(frame, NSamples, 2);
        % apply windowing to suprress sidelobes
        frame = frame .* windoppler3D .* winrange3D;
        % apply 3D FFT
        frameFFT = fftn(frame);
        % Range Doppler plot
        RD = fftshift(20*log10(squeeze(mean(abs(frameFFT),3))),1);
        frames(i, :, :) = RD(:, :);
        % Remove constant part by subtracting frame 1
        if i == 1
            frame1 = RD;
            fprintf("frame dimension %s\n", num2str(size(RD)));
        else
            difference = (RD - frame1 > 10);
            RD = RD .* difference;
        end
        if i == Nframes/2
            midframe = RD;
        end
        frames_nozero(i, :, :) = RD(:, :);
    end
    
    % Scale frames to match max range/doppler
    frames_scaled = frames(:, dop_bin_first:dop_bin_last, 1:length(ran_bin_mod));
    frames_nozero_scaled = frames_nozero(:, dop_bin_first:dop_bin_last, 1:length(ran_bin_mod));
    % fprintf("frames_scaled dim: %s, frames_nozero_scaled dim: %s\n", num2str(size(frames_scaled)), num2str(size(frames_nozero_scaled)));
    
    fprintf("Preprocessing data...\n");
    % Filter first
    frames_filtered = filterFrames(frames_nozero_scaled);
    
    % Convolution downsampling frames
    conv_mask = ones(9) / 9;
    frames_conv = zeros(Nframes, NDoppler, NRange);
    for i = 1:Nframes
        frames_conv(i, :, :) = conv2(squeeze(frames_filtered(i, :, :)), conv_mask, 'same');
    end
    
    fprintf("Finding cluster...\n");
    % Find cluster number
    max_k = 5;
    cutoff = 0.88; % tune this parameter
    [clusters, total_cluster] = findCluster(frames_conv, max_k, cutoff);
    % display(clusters);
    % display(total_cluster);
    fprintf("\nNumber of people detected: %d\n", total_cluster);
end

%% Functions
% ==============================================
% Plot one frame
% ==============================================
% ------------------------------------------------------------------
% Description: Plot one frame, dynamic scaling
%
% Method: Display the frame figure with the given axis
% Input:  Frame (2D), xaxis, yaxis
% ------------------------------------------------------------------
function plotOneFrame(frame, xaxis, yaxis, idx)
    figure;
    imagesc(xaxis, yaxis, frame);
    axis image;
    colorbar;
    xlabel('Range (m)');
    ylabel('Doppler Velocity (m/s)');
    plot_title = sprintf('Frame %d', idx);
    title(plot_title);
end

% ==============================================
% Plot multiple frames animation
% ==============================================
% ------------------------------------------------------------------
% Description: Plot all frames animation, dynamic scaling
%
% Method: Display figures sequentially with the given axis
% Input:  Frames (3D), xaxis, yaxis
% ------------------------------------------------------------------
function plotFrames(frames, xaxis, yaxis)
    figure(1);
    for i = 1:size(frames, 1)
        imagesc(xaxis, yaxis, squeeze(frames(i, :, :)));
        axis image;
        colorbar;
        xlabel('Range (m)');
        ylabel('Doppler Velocity (m/s)');
        pause(0.1);
    end
end

% ==============================================
% Frame filtering (all frames individually)
% ==============================================
% ------------------------------------------------------------------
% Description: Filter frames by their respective mean/median power
%
% Method: Extract mean and median power values for each frame and 
%         only keep the data points within certain threshold (set 
%         all else to zeros)
% Input:  Original frames (Nframes, NChirps, NSamples)
% Output: Filtered frames (Nframes, NChirps, NSamples)
% ------------------------------------------------------------------
function [frames_f] = filterFrames(frames)
    N_frame = size(frames, 1);
    N_chirp = size(frames, 2);
    N_sample = size(frames, 3);
    frames_f = zeros(N_frame, N_chirp, N_sample);
    for i = 1:size(frames,1)
        % compute average and max including zeros
        nonzero_mask = frames(i, :, :);
        average_power = mean(nonzero_mask, 'all');
        maximum_power = max(nonzero_mask, [], 'all');
        % compute medium excluding zeros
        nonzero_mask(nonzero_mask==0) = NaN;
        median_power = median(nonzero_mask, 'all', 'omitmissing');
        % Change threshold here
        lower_threshold = (average_power + median_power) / 2;
        upper_threshold = (median_power + maximum_power) / 12;
        % fprintf("frame %d, max = %f, avg = %f, med = %f, lower = %f, upper = %f\n", i, maximum_power, average_power, median_power, lower_threshold, upper_threshold);
        % filtering_mask = frames(i, :, :) >= lower_threshold & frames(i, :, :) <= upper_threshold;
        filtering_mask = frames(i, :, :) >= upper_threshold;
        frames_f(i, :, :) = frames(i, :, :) .* filtering_mask;
    end
end

% ==============================================
% Find aggregate cluster number from all frames
% ==============================================
% ------------------------------------------------------------------
% Description: Vote among all frames to determine overall k
% 
% Method: Count the most frequent k values in all frames
% Input:  All frames data (3D), Maximum k to test, Cutoff
% Output: Aggregated k value
% ------------------------------------------------------------------
function [Ks, K] = findCluster(X, varargin)
    if nargin > 1, max_k = cell2mat(varargin(1)); else max_k = 3; end
    if nargin > 2, cutoff = cell2mat(varargin(2)); else cutoff = 0.95; end

    N_frame = size(X,1);
    Ks = zeros(1, N_frame);
    for i = 1:N_frame
        [one_k, ~] = kmeansOptk(squeeze(X(i, :, :)), max_k, cutoff);
        Ks(i) = one_k;
    end
    K = mode(Ks);
end

% =================================================
% K-means clustering with Elbow Method (per frame)
% =================================================
% ------------------------------------------------------------------
% Description: Find the optimal cluster number for each frame
%
% Method: Compute K-means Within-Cluster-Sum of Squared Errors (WSS)
%         through iteration on different K values. Find the optimal
%         k that pass the cutoff of variance.
% Input:  Single frame data (2D), Maximum k to test, Cutoff
% Output: Optimal cluster number k and the optimal distance SUMD
% Reference: Sebastien De Landtsheer (2023). kmeans_opt (https://www.mathworks.com/matlabcentral/fileexchange/65823-kmeans_opt), 
%            MATLAB Central File Exchange. Retrieved June 3, 2023.
% ------------------------------------------------------------------
function [K, SUMD] = kmeansOptk(X, varargin)
    if nargin > 1, MaxK = cell2mat(varargin(1)); else MaxK = 5; end
    if nargin > 2, Cutoff = cell2mat(varargin(2)); else Cutoff = 0.95; end
    % fprintf("MaxK: %d, Cutoff: %d\n", MaxK, Cutoff)

    D = zeros(MaxK);
    % For each k value, get the minimum WSS
    for k = 1:MaxK
        [~, ~, dist] = kmeans(X, k);
        tmp = sum(dist);
        % Repeat kmeans algorithm for randomization purpose
        for cc=2:3
            [~, ~, dist] = kmeans(X, k);
            tmp = min(sum(dist), tmp);
        end
        D(k) = tmp;
    end
    
    % Compute diff between each WSS
    Var = D(1:end-1) - D(2:end);
    PC = cumsum(Var) / (D(1) - D(end));
    
    r = find(PC > Cutoff);
    K = r(1);
    [~, ~, SUMD] = kmeans(X, K);
end


%% Reset all
clear all;
clc;
close all;

%% Import radarcube
data = load("C:\Users\Owner\Downloads\dataCube-20230606T052836Z-001\dataCube\radarCube_0531_8.mat"); % cells are frames, each cell is (NChirps, NReceivers, NSamples)

%% Bin ranges
radarCube = data.radarCube;

numRangeBins = radarCube.rfParams.numRangeBins;
numChirps = radarCube.dim.numChirps;
rangeResolutionsInMeters = radarCube.rfParams.rangeResolutionsInMeters;
dopplerResolutionMps = radarCube.rfParams.dopplerResolutionMps;

ran_bin = (1:numRangeBins) * rangeResolutionsInMeters;
% dop_bin = (1:d_bin) * dopplerResolutionMps;
tot_dop_bin = ((numChirps/2) - (1:numChirps)) * dopplerResolutionMps;

%% Apply windowing and FFT processing
Nframes = length(radarCube.data);
[NChirps, NReceivers, NSamples] = size(radarCube.data{1});
fprintf("NChirps: %d, NReceivers: %d, NSamples: %d\n", NChirps, NReceivers, NSamples);

%Process each frame
winrange = hann(NSamples);
winrange3D = permute(repmat(winrange, [1, NChirps, NReceivers]), [2, 1, 3]);
windoppler = hann(NChirps);
windoppler3D = repmat(windoppler, [1, NSamples, NReceivers]);

average_ppl_detected = 0;

% Individual frames to save
frame1 = zeros(NChirps, NSamples);
midframe = zeros(NChirps, NSamples);
% Initialize cell arrays to store peak locations for each frame
range_values_per_frame = cell(Nframes, 1);
doppler_values_per_frame = cell(Nframes, 1);
figure(1);
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
    % Capping at 4m Max Range [Customizable input]
    range_indices = ran_bin < 4;
    ran_bin_mod = ran_bin(range_indices);
    % Capping Doppler between -2, 2:
    doppler_indices = (tot_dop_bin >= -2) & (tot_dop_bin <= 2);
    tot_dop_bin_mod = tot_dop_bin(doppler_indices);
    % Filtered RD:
    RD_filtered = RD(95:161, 1:length(ran_bin_mod));
    % imagesc(ran_bin_mod, tot_dop_bin_mod, RD_filtered);
    % colorbar;
    % axis image;
    % clim([70, 125]); % adjust min and max intensity to match your data
    % xlabel('Range (m)');
    % ylabel('Doppler Velocity (m/s)');
    % pause(0.1); 

    peak_threshold = 0.95 * max(RD_filtered(:));
    % Find the peaks in the filtered data
    [peak_values, peak_indices] = findpeaks(RD_filtered(:));
    peak_indices = peak_indices(peak_values > peak_threshold);
    peak_values = peak_values(peak_values > peak_threshold);
    
    % Convert peak indices to subscripts
    [x, y] = ind2sub(size(RD_filtered), peak_indices);
    range_values = ran_bin_mod(y);
    doppler_values = tot_dop_bin_mod(x);
    % Store range and doppler values for the current frame
    range_values_per_frame{i} = range_values;
    doppler_values_per_frame{i} = doppler_values;


    [range_values doppler_values] = cluster(range_values, doppler_values);
    test = numel(range_values)
    average_ppl_detected = average_ppl_detected + numel(range_values);
    
    % Plot the doppler-range plot with detected peaks
    imagesc(ran_bin_mod, tot_dop_bin_mod, RD_filtered);
    hold on;
    plot(range_values, doppler_values, 'rx', 'MarkerSize', 10);
    colorbar;
    xlabel('Range (m)');
    ylabel('Doppler Velocity (m/s)');
    pause(0.1); 
end

num_ppl = average_ppl_detected / Nframes
num_ppl = round(num_ppl)


%% Clustering
function [range_values, doppler_values] = cluster(range_values, doppler_values)
dummy_range_values = range_values;
dummy_doppler_values = doppler_values;

for i = 1:(numel(range_values) - 1.5)
    for j = (i + 1):numel(range_values)
       pix_1 = [range_values(i), doppler_values(i)];
       pix_2 = [range_values(j), doppler_values(j)];
       distance = sqrt((pix_1(1)-pix_2(1))^2 + (pix_1(2)-pix_2(2))^2);

       if (distance < 0.15) % 0.3 is the best value so far
          dummy_doppler_values(j) = 0;
          dummy_range_values(j) = 0;
       end

    end
end

doppler_values = nonzeros(dummy_doppler_values');
range_values = nonzeros(dummy_range_values');

return
end

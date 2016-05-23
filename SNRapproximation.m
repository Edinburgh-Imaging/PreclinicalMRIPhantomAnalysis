function [level, meansignal, sdnoise, SNR] = SNRapproximation(vol)
% Companion function of the PhantomAnalysisGUI. It measures the
% signal-to-noise ratio (SNR) of a greyscale 3D matrix using approximate
% distributions of the signal and noise. This assumes that the two 
% distributions are distinguishable (high contrast) and that the noise is
% described by a Rician distribution with a higher peak than the signal's
% peak.
%
% >> [level, meansignal, sdnoise, SNR] = SNRapproximation(vol)
%
% Variable Dictionary:
% --------------------
% vol           input     The greyscale 3D matrix.
% level         output    The intensity level splitting the distributions.
% meansignal    output    The mean intensity of the 'high' signal.
% sdnoise       output    The SD of the intensity of the 'low' signal.
% SNR           output    The SNR.
%
% Last Modified: 25 January 2016
% Copyright (c) 2016, Xenios Milidonis

% Get the histogram of the whole matrix with 500 bins.
[pixelcounts, greylevels] = imhist(vol(:), 500);

% The tallest peak should belong to the background noise. If the peak is at
% positions 0 or 500, exclude them from the histogram by setting them to 0.
[~, maxindex] = max(pixelcounts);
if maxindex == 1                
    pixelcounts(1) = 0;
    [~, maxindex] = max(pixelcounts);
end
if maxindex == 500                
    pixelcounts(500) = 0;
    [~, maxindex] = max(pixelcounts);
end

% The intensity threshold to split the distributions should be about 5 
% times the position of the peak.
level = 5 * greylevels(maxindex);

% Get the two distributions above and below the threshold.
signal = vol(vol > level);
noise = vol(vol <= level);

% Estimate the mean signal and the SD of the noise.
meansignal = mean(signal);
sdnoise = std(noise);

% Estimate the SNR using a factor to account for the Rician noise 
% distribution in MRI magnitude images.
SNR = 0.655 * meansignal / sdnoise;


function vol = stacktomatrix(image)
% It transforms a TIFF stack of images into an intensity-normalised 3D
% matrix of type double.
%
% >> vol = stacktomatrix(image)
%
% Variable Dictionary:
% --------------------
% image    input     The path of the TIFF image sequence/stack.
% vol      output    The created 3D matrix.
%
% Last Modified: 16 March 2016
% Copyright (c) 2016, Xenios Milidonis

% Get the slice number.
info = imfinfo(image);
numofslices = numel(info);

% Create a 3D matrix version of the input image sequence.
for i = 1:numofslices                                              
    slice = imread(image, 'Index', i, 'Info', info);
    slice = slice(:, :, 1); 
    
    slice = uint16(slice);
    slice = im2double(slice);

    if i == 1                                       
        vol = slice;
    else
        vol = cat(3, vol, slice);
    end
end

% Normalise the image intensities from 0 to 1.
vol = vol - min(vol(:));
vol = vol / max(vol(:));


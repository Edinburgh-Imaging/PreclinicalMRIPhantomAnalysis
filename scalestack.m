function scaledseq = scalestack(vol, ni, nf, scalefactor, method)
% Companion function of the PhantomAnalysisGUI. Specified slices through
% the input 3D matrix are scaled according to the input scale factor and
% the specified resizing method. A new 3D matrix is created.
%
% >> scaledseq = scalestack(vol, ni, nf, scalefactor, method)
%
% Variable Dictionary:
% --------------------
% vol           input     The original 3D matrix.
% ni            input     Start from this slice.
% nf            input     Finish at this slice.
% scalefactor   input     Scale slices by this factor.
% method        input     The method to resize the 3D matrix: 'nearest',
%                         'bilinear' or 'bicubic'.
% scaledseq     output    The scaled 3D matrix.
%
% Last Modified: 31 January 2016
% Copyright (c) 2016, Xenios Milidonis

for i = ni:nf                             
    % Resize using the default bicubic interpolation.
    scaledset = imresize(vol(:, :, i), scalefactor, method);  

    if i == ni
        scaledseq = scaledset;
    else
        scaledseq = cat(3, scaledseq, scaledset);
    end
end
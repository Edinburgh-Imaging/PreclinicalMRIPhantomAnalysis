function bwimage = kmeanssegmentation(image, k, l)
% Companion function of the PhantomAnalysisGUI. It segments a grayscale 
% matrix using k-means clustering, with the l brightest clusters classified
% as the region of interest out of a total of k clusters.
%
% >> bwimage = kmeanssegmentation(image, k, l)
%
% Variable Dictionary:
% --------------------
% image        input    The greyscale image or volume (any dimensions).
% k            input    Number of total clusters required.
% l            input    Number of clusters that should be classified as  
%                       the region of interest. The 'brightest' classes are 
%                       chosen and the rest are classified as background.
%                       l must be larger than 0 and smaller than k.
% bwimage      output   The segmented binary image or volume.
%
% Last Modified: 18 February 2016
% Copyright (c) 2016, Xenios Milidonis

% Check input variables.
if (l == 0) || (l >= k) 
    error('l must be larger than 0 and smaller than k.');
end

% k-means requires double precision.
image = im2double(image);

% Cluster. Pixels in the clustered image will have values from 1 to k.
clustered = kmeans(image(:), k, 'EmptyAction', 'singleton');

% The output is a [n x m, 1] column vector so it must be reshaped.
clustered = reshape(clustered, size(image));

% k-means initializes at random locations in the image thus clusters will 
% not be sorted according to their intensities by default. Sort by first 
% creating a vector with the intensity of the first pixel in each cluster.
clusterintensity = zeros(k, 1);
for i = 1:k
    clusterintensity(i) = image(find(clustered == i, 1));
end 

% Then, find the index of each cluster in order.
clusteridx = zeros(k, 1);
for i = 1:k
    % Find the minimum intensity; set this as the first cluster.
    clusteridx(find(clusterintensity == min(clusterintensity))) = i; 
    
    % Set the intensity of the first cluster in the vector as NaN in order
    % to find the next minimum intensity.
    clusterintensity(clusterintensity == min(clusterintensity)) = NaN;
end

% Finally, rearrange the clustered image.
bwimage = clusteridx(clustered);

% Make the l brightest clusters white and the rest black.
bwimage(bwimage <= (k - l)) = 0; 
bwimage(bwimage > (k - l)) = 1;


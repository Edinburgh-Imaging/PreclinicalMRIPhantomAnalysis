function varargout = PhantomAnalysisGUI(varargin)
% PHANTOMANALYSISGUI MATLAB code for PhantomAnalysisGUI.fig
% This is a GUI used for the automated analysis of MRI images of a small 
% structural phantom for the assessment of geometric accuracy in 
% preclinical MRI scanners.
%
% Use the GUI:
% Please see the README text file.
%
% Edit the GUI:
% >> guide
% 
% Last Modified: 11 August 2016
% Copyright (c) 2016, Xenios Milidonis

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PhantomAnalysisGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @PhantomAnalysisGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% _________________________________________________________________________
% --- Executes just before PhantomAnalysisGUI is made visible.
function PhantomAnalysisGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.

global imagefolder;

% Choose default command line output for PhantomAnalysisGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Set the current directory as the starting path to look for the images.
imagefolder = pwd;

% Set the current date and time as the first output to the GUI's log window.
set(handles.listbox_log, 'String', datestr(now));

% Move the GUI to the centre of the screen.
movegui(hObject, 'center');

% Hide the axes values from all axes objects.
set(handles.axes_inputimage, 'XTick', []);
set(handles.axes_inputimage, 'YTick', []);
set(handles.axes_outputimage, 'XTick', []);
set(handles.axes_outputimage, 'YTick', []);

% Set some default parameters.
set(handles.radiobutton_2clusters, 'Value', 0);
set(handles.radiobutton_3clusters, 'Value', 1);
set(handles.radiobutton_kmeansdim, 'Value', 0);
set(handles.radiobutton_cannydim, 'Value', 1);
set(handles.radiobutton_kmeansvol, 'Value', 0);
set(handles.radiobutton_cannyvol, 'Value', 1);
set(handles.popup_show, 'Value', 1);

% Global variables must be cleared, otherwise they remain in workspace.
clear global imagename;
clear global numofslices;
clear global centralslice;
clear global vol;
clear global volaligned;
clear global voladjusted;
clear global maxvaluevol;
clear global ipsimask;
clear global contramask;
clear global analyseddimmasks;
clear global outlinedimmasks;
clear global analysedvolmasks;
clear global outlinevolmasks;
clear global locationsArray;

% UIWAIT makes PhantomAnalysisGUI wait for user response (see UIRESUME)
% uiwait(handles.phantomanalysisfigure);

% --- Outputs from this function are returned to the command line.
function varargout = PhantomAnalysisGUI_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;


% _________________________________________________________________________
% --- Executes on button press in pushbutton_loaddata.
function pushbutton_loaddata_Callback(hObject, eventdata, handles)

global imagefolder;
global imagename;
global numofslices;
global centralslice;
global vol;
global volaligned;
global voladjusted;
global maxvaluevol;

warning('off', 'all')

% Open UI to select input image stack and updade imagefolder.
[imagename, imagefolder] = uigetfile(fullfile(imagefolder, '*.tif'),... 
    'Please select a TIFF image sequence');
imagepath = fullfile(imagefolder, imagename);

% Get the number of slices in the image.
info = imfinfo(imagepath);
numofslices = numel(info);

% Show loaded image path in text_imagename.
set(handles.text_imagename, 'String', imagepath);

% -----------------ADJUST THE DATA AND ESTIMATE PARAMETERS-----------------
% Call stacktomatrix.m to transform the TIFF stack into a 3D matrix. This
% will be passed to slider_inputslices and uipanel_plot to plot the 
% histogram of the original data.
vol = stacktomatrix(imagepath);
     
% Call identifyslice.m to find the first, central and last slice.
[firstslice, centralslice, lastslice] = identifyslice(vol);
central = vol(:, :, centralslice);

% Create an aligned volume to be passed for calculations and another
% contrast-adjusted and aligned volume to be passed for visualisation in
% the GUI's 'Output Data' figure.
volaligned  = vol;
voladjusted = zeros(size(vol));

% Adjust the contrast of voladjusted.
for i = 1:numofslices     
    voladjusted(:, :, i) = imadjust(vol(:, :, i));
end

% Call phantomorientation.m to get the angle to align the images.
alignangle = phantomorientation(central, 'b');

% Rotate using the align angle.
volaligned  = imrotate(volaligned,  alignangle, 'bicubic', 'crop');  
voladjusted = imrotate(voladjusted, alignangle, 'bicubic', 'crop');

% Show a warning message in the GUI's log window stating this angle...
oldmsgs = cellstr(get(handles.listbox_log, 'String'));
newmsg = sprintf('Phantom misaligned by %4.1f degrees.', alignangle);
set(handles.listbox_log, 'String', [oldmsgs; newmsg] );
%     ... and make its scroll bar move at the bottom.
logsize = size(get(handles.listbox_log, 'String'), 1);
set(handles.listbox_log, 'Value', logsize);

% Show the central slice of the stack with contrast adjustment.
axes(handles.axes_inputimage);
imshow(imadjust(vol(:, :, centralslice)));
set(handles.axes_inputimage, 'XTick', []);
set(handles.axes_inputimage, 'YTick', []);

% Call kmeanssegmentation.m to get the binarised central slice using 
% k-means clustering. Two white clusters out of a total of 3 is the best 
% option for phantom images with no uniform illumination.
bwcentral = kmeanssegmentation(central, 3, 2);

% If multiple objects are found, try to attach them.
bwcentral = attachobjects(bwcentral);

% Call centroid.m to get the position of the centroid and round it up.
[centroidr, centroidc] = centroid(bwcentral);
centroidr = round(centroidr);
centroidc = round(centroidc);

% Get the imaging plane.
horvector = bwcentral(centroidr, :);
vervector = bwcentral(:, centroidc);
horpixelcount = sum(horvector);
verpixelcount = sum(vervector);
if horpixelcount < (0.85 * verpixelcount)
    plane = 'Sagittal'; % Or Axial90. Rotate.
elseif horpixelcount > (1.15 * verpixelcount)
    plane = 'Axial';    % Or Sagittal90. Do not rotate.
else 
    plane = 'Coronal';  % Or Coronal90. Do not rotate.
end

% Rotate images by 90 degrees if required.
if strcmp(plane, 'Sagittal')
    volaligned  = imrotate(volaligned,  -90);
    voladjusted = imrotate(voladjusted, -90);
end 

% Call SNRapproximation.m to get an approximation of the background noise 
% and signal in the data for estimating the sigma for Canny edge detector.
% The size of the images is also considered.
[~, meansignal, sdnoise, ~] = SNRapproximation(vol);
sigma = size(vol, 1) / (30 * sqrt(meansignal / sdnoise));

% If SNRapproximation.m fails to calculate a sigma (NaN), use a starting
% value.
if isnan(sigma)
    sigma = 1.5;
end

% ----------------------SET DEFAULT PARAMETERS ON GUI----------------------
% Show the central slice number in text_inputslicenumber.
set(handles.text_inputslicenumber, 'String',...
    sprintf('%d/%d', centralslice, numofslices));

% Set image slider's initial position, min, max and step.
set(handles.slider_inputslices, 'Value', centralslice);
set(handles.slider_inputslices, 'Min', 1);
set(handles.slider_inputslices, 'Max', numofslices);
if numofslices == 1
    set(handles.slider_inputslices, 'SliderStep', [1, 1]);
else
    set(handles.slider_inputslices, 'SliderStep',...
        [1 / (numofslices - 1), 1 / (numofslices - 1)]);
end

% Show image's histogram (default is stack histogram).
set(handles.radiobutton_stack, 'Value', 1);
set(handles.radiobutton_slice, 'Value', 0);
maxvaluevol = max(vol(:));
axes(handles.axes_histogram);
imhist(vol(:), 500);
axis([0 maxvaluevol 0 inf])
set(handles.axes_histogram, 'FontSize', 7);

% Set the image resolution, if found in the header.
if isfield(info(1), 'XResolution') && ~isempty(info(1).XResolution)
    xres = info(1).XResolution;
    hpixsize = 1/xres;
    set(handles.text_xres, 'String', hpixsize);
end
if isfield(info(1), 'YResolution') && ~isempty(info(1).YResolution)
    yres = info(1).YResolution;
    vpixsize = 1/yres;
    set(handles.text_yres, 'String', vpixsize);
end

% Set inputs for analysis.
if strcmp(plane, 'Axial')
    set(handles.popup_plane, 'Value', 1);
elseif strcmp(plane, 'Coronal')
    set(handles.popup_plane, 'Value', 2);
else
    set(handles.popup_plane, 'Value', 3);
end
set(handles.text_firstslice, 'String', firstslice);
set(handles.text_lastslice, 'String', lastslice);
set(handles.text_scalefactor, 'String', 5);
set(handles.text_edgethreshold, 'String', 2.5);
set(handles.text_edgesigma, 'String', sprintf('%.2f', sigma));


% _________________________________________________________________________
% --- Executes on slider movement.
function slider_inputslices_Callback(hObject, eventdata, handles)

global vol;
global numofslices;

% Get the current slice number from the slider.
sliderslice = round(get(hObject, 'Value'));

% Show the contrast-adjusted image. The voladjusted is not used as it is 
% rotated; the original orientation must be shown.
axes(handles.axes_inputimage);
imshow(imadjust(vol(:, :, sliderslice)));
set(handles.axes_inputimage, 'XTick', []);
set(handles.axes_inputimage, 'YTick', []);

% Plot the histogram of individual slices of the original stack.
histogrambutton = get(handles.radiobutton_slice, 'Value');
axes(handles.axes_histogram);
if histogrambutton == 1
    maxvalueslice = max(max(vol(:, :, sliderslice)));
    imhist(vol(:, :, sliderslice), 500);
    axis([0 maxvalueslice 0 inf])
    set(handles.axes_histogram, 'FontSize', 7);
end

% Show the slice number in text_inputslicenumber
set(handles.text_inputslicenumber, 'String',...
    sprintf('%d/%d', sliderslice, numofslices));

% --- Executes during object creation, after setting all properties.
function slider_inputslices_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% _________________________________________________________________________
% --- Executes when selected object is changed in uipanel_plot.
function uipanel_plot_SelectionChangeFcn(hObject, eventdata, handles)

global vol;
global maxvaluevol;

% Get the tag of the selected button.
histogrambutton = get(eventdata.NewValue, 'Tag');

% If the stack button is selected then plot its histogram, otherwise rely 
% on slider_inputslices to plot individual slice histograms.
axes(handles.axes_histogram);
if strcmp(histogrambutton, 'radiobutton_stack')
    imhist(vol(:), 500);
    axis([0 maxvaluevol 0 inf])
    set(handles.axes_histogram, 'FontSize', 7);
end

% INDIVIDUAL BUTTONS OF A BUTTON GROUP MUST NOT BE CODED
% --- Executes during object creation, after setting all properties.
function radiobutton_stack_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function radiobutton_slice_CreateFcn(hObject, eventdata, handles)


% _________________________________________________________________________
% --- Executes on button press in radiobutton_2clusters.
function radiobutton_2clusters_Callback(hObject, eventdata, handles)

% Release radiobutton_3clusters.
set(handles.radiobutton_3clusters, 'Value', 0);


% _________________________________________________________________________
% --- Executes on button press in radiobutton_3clusters.
function radiobutton_3clusters_Callback(hObject, eventdata, handles)

% Release radiobutton_2clusters.
set(handles.radiobutton_2clusters, 'Value', 0);


% _________________________________________________________________________
% --- Executes on button press in radiobutton_kmeansdim.
function radiobutton_kmeansdim_Callback(hObject, eventdata, handles)

% Release radiobutton_cannydim.
set(handles.radiobutton_cannydim, 'Value', 0);


% _________________________________________________________________________
% --- Executes on button press in radiobutton_cannydim.
function radiobutton_cannydim_Callback(hObject, eventdata, handles)

% Release radiobutton_kmeansdim.
set(handles.radiobutton_kmeansdim, 'Value', 0);


% _________________________________________________________________________
% --- Executes on button press in radiobutton_kmeansvol.
function radiobutton_kmeansvol_Callback(hObject, eventdata, handles)

% Release radiobutton_cannyvol.
set(handles.radiobutton_cannyvol, 'Value', 0);


% _________________________________________________________________________
% --- Executes on button press in radiobutton_cannyvol.
function radiobutton_cannyvol_Callback(hObject, eventdata, handles)

% Release radiobutton_kmeansvol.
set(handles.radiobutton_kmeansvol, 'Value', 0);


% ______________________________ ANALYSIS _________________________________
% _________________________________________________________________________
% --- Executes on button press in checkbox_saveimages.
function checkbox_saveimages_Callback(hObject, eventdata, handles)


% _________________________________________________________________________
% --- Executes on button press in pushbutton_dimensions.
function pushbutton_dimensions_Callback(hObject, eventdata, handles)

set(handles.togglebutton_action, 'Value', 1);

global imagename;
global imagefolder;
global centralslice;
global volaligned;
global voladjusted;
global analyseddimmasks;
global locationsArray;

% Get all required parameter values. If a value is not given, make the
% corresponding GUI field red.
planevalue = get(handles.popup_plane, 'Value');
planestring = get(handles.popup_plane, 'String');
plane = planestring{planevalue};

if isempty(get(handles.text_xres, 'String'))
    set(handles.text_xres, 'BackgroundColor', 'r');
else
    set(handles.text_xres, 'BackgroundColor', 'w');
    hpixsize = str2double(get(handles.text_xres, 'String'));
end
if isempty(get(handles.text_yres, 'String'))
    set(handles.text_yres, 'BackgroundColor', 'r');
else
    set(handles.text_yres, 'BackgroundColor', 'w');
    vpixsize = str2double(get(handles.text_yres, 'String'));
end
if isempty(get(handles.text_firstslice, 'String'))
    set(handles.text_firstslice, 'BackgroundColor', 'r');
else
    set(handles.text_firstslice, 'BackgroundColor', 'w');    
    ni = str2double(get(handles.text_firstslice, 'String'));
end
if isempty(get(handles.text_lastslice, 'String'))
    set(handles.text_lastslice, 'BackgroundColor', 'r');
else
    set(handles.text_lastslice, 'BackgroundColor', 'w');    
    nf = str2double(get(handles.text_lastslice, 'String'));
end
if isempty(get(handles.text_scalefactor, 'String'))
    set(handles.text_scalefactor, 'BackgroundColor', 'r');
else
    set(handles.text_scalefactor, 'BackgroundColor', 'w');    
    scalefactor = str2double(get(handles.text_scalefactor, 'String'));
end
if isempty(get(handles.text_edgethreshold, 'String'))
    set(handles.text_edgethreshold, 'BackgroundColor', 'r');
else
    set(handles.text_edgethreshold, 'BackgroundColor', 'w');    
    thresholdf = str2double(get(handles.text_edgethreshold, 'String'));
end
if isempty(get(handles.text_edgesigma, 'String'))
    set(handles.text_edgesigma, 'BackgroundColor', 'r');
else
    set(handles.text_edgesigma, 'BackgroundColor', 'w');    
    sigma = str2double(get(handles.text_edgesigma, 'String'));
end
set(handles.text_slicethickness, 'BackgroundColor', 'w');

% Check if user wants to save images. This equals 1 if yes, 0 if not.
saveimages = get(handles.checkbox_saveimages, 'Value');

% Calculate the dimensions of the phantom with the method asked (k-means or
% Canny edge detection) and get the boundary-only masks.
methodbutton = get(handles.radiobutton_kmeansdim, 'Value');
if methodbutton == 1       % Use k-means.
    clustersbutton = get(handles.radiobutton_2clusters, 'Value');
    if clustersbutton == 1 % Use k = 2.
        method = 'kmeans2';
    else                   % Use k = 3.
        method = 'kmeans3';
    end
else                       % Use Canny.
    method = 'Canny';
end

% Call phantomdimensions.m to estimate the phantom dimensions.
[analyseddimmasks, locationsArray, xdimension, ydimension] =...
    phantomdimensions(imagename, imagefolder, volaligned, plane, hpixsize,...
    vpixsize, centralslice, ni, nf, scalefactor, thresholdf, sigma,...
    saveimages, method);

% Set the calculated dimensions in the respective edit boxes.
set(handles.text_xdimension, 'String', xdimension);
set(handles.text_ydimension, 'String', ydimension);

% Show the first slice of the analysed image.
showvalue = get(handles.popup_show, 'Value');
axes(handles.axes_outputimage);
if showvalue == 2  % masks only
    imshow(analyseddimmasks(:, :, 1));
else               % overlaid
    imshow(voladjusted(:, :, ni));
    red = cat(3, ones(size(voladjusted(:, :, ni))),...
        zeros(size(voladjusted(:, :, ni))),...
        zeros(size(voladjusted(:, :, ni))));
    hold on
    h = imshow(red);
    hold off
    set(h, 'AlphaData', analyseddimmasks(:, :, 1));
    for i = 1:16
        hold on
        plot(locationsArray(i, 2, 1), locationsArray(i, 1, 1),...
            'xg', 'MarkerSize', 8);
        hold off
    end
end
set(handles.axes_outputimage, 'XTick', []);
set(handles.axes_outputimage, 'YTick', []);

% Show the slice number in text_outputslicenumber.
set(handles.text_outputslicenumber, 'String', ni);

% Set analysed image slider's initial position, min, max and step.
set(handles.slider_outputslices, 'Value', ni);
set(handles.slider_outputslices, 'Min', ni);
set(handles.slider_outputslices, 'Max', nf);
if ni == nf
    set(handles.slider_outputslices, 'SliderStep', [1, 1]);
else
    set(handles.slider_outputslices, 'SliderStep',...
        [1 / (nf - ni), 1 / (nf - ni)]);
end

% Output to GUI's log window.
oldmsgs = cellstr(get(handles.listbox_log, 'String'));
newmsgx = sprintf('%s: %8.3f mm horizontally', method, xdimension);
newmsgy = sprintf('%s: %8.3f mm vertically', method, ydimension);
newmsgx = ['<HTML><FONT color = "blue"><b>', newmsgx]; % make it bold blue
newmsgy = ['<HTML><FONT color = "blue"><b>', newmsgy]; % make it bold blue
set(handles.listbox_log, 'String', [oldmsgs; newmsgx; newmsgy]);
logsize = size(get(handles.listbox_log, 'String'), 1);
set(handles.listbox_log, 'Value', logsize);


% _________________________________________________________________________
% --- Executes on button press in pushbutton_volume.
function pushbutton_volume_Callback(hObject, eventdata, handles)

set(handles.togglebutton_action, 'Value', 0);

global imagename;
global imagefolder;
global centralslice;
global volaligned;
global voladjusted;
global analysedvolmasks;
global outlinevolmasks;

% Get all required parameter values. If a value is not given, make the
% corresponding GUI field red.
planevalue = get(handles.popup_plane, 'Value');
planestring = get(handles.popup_plane, 'String');
plane = planestring{planevalue};

if isempty(get(handles.text_xres, 'String'))
    set(handles.text_xres, 'BackgroundColor', 'r');
else
    set(handles.text_xres, 'BackgroundColor', 'w');
    hpixsize = str2double(get(handles.text_xres, 'String'));
end
if isempty(get(handles.text_yres, 'String'))
    set(handles.text_yres, 'BackgroundColor', 'r');
else
    set(handles.text_yres, 'BackgroundColor', 'w');
    vpixsize = str2double(get(handles.text_yres, 'String'));
end
if isempty(get(handles.text_slicethickness, 'String'))
    set(handles.text_slicethickness, 'BackgroundColor', 'r');
else
    set(handles.text_slicethickness, 'BackgroundColor', 'w');
    slicethickness = str2double(get(handles.text_slicethickness, 'String'));
end
if isempty(get(handles.text_firstslice, 'String'))
    set(handles.text_firstslice, 'BackgroundColor', 'r');
else
    set(handles.text_firstslice, 'BackgroundColor', 'w');    
    ni = str2double(get(handles.text_firstslice, 'String'));
end
if isempty(get(handles.text_lastslice, 'String'))
    set(handles.text_lastslice, 'BackgroundColor', 'r');
else
    set(handles.text_lastslice, 'BackgroundColor', 'w');    
    nf = str2double(get(handles.text_lastslice, 'String'));
end
if isempty(get(handles.text_edgethreshold, 'String'))
    set(handles.text_edgethreshold, 'BackgroundColor', 'r');
else
    set(handles.text_edgethreshold, 'BackgroundColor', 'w');    
    thresholdf = str2double(get(handles.text_edgethreshold, 'String'));
end
if isempty(get(handles.text_edgesigma, 'String'))
    set(handles.text_edgesigma, 'BackgroundColor', 'r');
else
    set(handles.text_edgesigma, 'BackgroundColor', 'w');    
    sigma = str2double(get(handles.text_edgesigma, 'String'));
end
set(handles.text_scalefactor, 'BackgroundColor', 'w');

% Check if user wants to save images. This equals 1 if yes, 0 if not.
saveimages = get(handles.checkbox_saveimages, 'Value');

% Calculate the volume of the cylindrical phantom compartment with the 
% method asked (k-means or Canny edge detection).
methodbutton = get(handles.radiobutton_kmeansvol, 'Value');
if methodbutton == 1       % Use k-means.
    clustersbutton = get(handles.radiobutton_2clusters, 'Value');
    if clustersbutton == 1 % Use k = 2.
        method = 'kmeans2';
    else                   % Use k = 3.
        method = 'kmeans3';
    end
else                       % Use Canny.
    method = 'Canny';
end

% Call phantomvolume.m to estimate the phantom volume.
[analysedvolmasks, volume] = phantomvolume...
    (imagename, imagefolder, volaligned, plane, hpixsize, vpixsize,...
    slicethickness, centralslice, ni, nf, thresholdf, sigma,...
    saveimages, method);

% Set the calculated volume in the respective edit box.
set(handles.text_volume, 'String', volume);

% Create a boundary-only version of the masks.
outlinevolmasks = bwperim(analysedvolmasks, 4);

% Show the first slice of the analysed image.
showvalue = get(handles.popup_show, 'Value');
axes(handles.axes_outputimage);
if showvalue == 2  % masks only
    imshow(analysedvolmasks(:, :, 1));
else               % overlaid
    imshow(voladjusted(:, :, ni));
    red = cat(3, ones(size(voladjusted(:, :, ni))),...
        zeros(size(voladjusted(:, :, ni))),...
        zeros(size(voladjusted(:, :, ni))));
    hold on
    h = imshow(red);
    hold off
    set(h, 'AlphaData', outlinevolmasks(:, :, 1));
end
set(handles.axes_outputimage, 'XTick', []);
set(handles.axes_outputimage, 'YTick', []);

% Show the slice number in text_outputslicenumber.
set(handles.text_outputslicenumber, 'String', ni);

% Set analysed image slider's initial position, min, max and step.
set(handles.slider_outputslices, 'Value', ni);
set(handles.slider_outputslices, 'Min', ni);
set(handles.slider_outputslices, 'Max', nf);
if ni == nf
    set(handles.slider_outputslices, 'SliderStep', [1, 1]);
else
    set(handles.slider_outputslices, 'SliderStep',...
        [1 / (nf - ni), 1 / (nf - ni)]);
end

% Output to GUI's log window.
oldmsgs = cellstr(get(handles.listbox_log, 'String'));
newmsg = sprintf('%s: %9.3f mm^3', method, volume);
newmsg = ['<HTML><FONT color = "blue"><b>', newmsg]; % make it bold blue
set(handles.listbox_log, 'String', [oldmsgs; newmsg] );
logsize = size(get(handles.listbox_log, 'String'), 1);
set(handles.listbox_log, 'Value', logsize);


% _________________________________________________________________________
% --- Executes on slider movement.
function slider_outputslices_Callback(hObject, eventdata, handles)

global voladjusted;
global analyseddimmasks;
global analysedvolmasks;
global outlinevolmasks;
global locationsArray;

% Get the number of the first analysed slice.
ni = str2double(get(handles.text_firstslice, 'String'));

% Get the current slice number from the slider.
sliderslice = round(get(hObject, 'Value'));

% Choose between analyseddim and analysedvol and show the analysed image.
actionvalue = get(handles.togglebutton_action, 'Value');
showvalue = get(handles.popup_show, 'Value');

axes(handles.axes_outputimage);
if (actionvalue == 1) && (showvalue == 2)      % masks only
    imshow(analyseddimmasks(:, :, sliderslice-ni+1));
elseif (actionvalue == 1) && (showvalue == 1)  % overlaid
    imshow(voladjusted(:, :, sliderslice));
    red = cat(3, ones(size(voladjusted(:, :, sliderslice))),...
        zeros(size(voladjusted(:, :, sliderslice))),...
        zeros(size(voladjusted(:, :, sliderslice))));
    hold on
    h = imshow(red);
    hold off
    set(h, 'AlphaData', analyseddimmasks(:, :, sliderslice-ni+1));
    for i = 1:16
        hold on
        plot(locationsArray(i, 2, sliderslice-ni+1),...
            locationsArray(i, 1, sliderslice-ni+1), 'xg', 'MarkerSize', 8);
        hold off
    end
elseif (actionvalue == 0) && (showvalue == 2)  % masks only
    imshow(analysedvolmasks(:, :, sliderslice-ni+1));
elseif (actionvalue == 0) && (showvalue == 1)  % overlaid
    imshow(voladjusted(:, :, sliderslice));
    red = cat(3, ones(size(voladjusted(:, :, sliderslice))),...
        zeros(size(voladjusted(:, :, sliderslice))),...
        zeros(size(voladjusted(:, :, sliderslice))));
    hold on
    h = imshow(red);
    hold off
    set(h, 'AlphaData', outlinevolmasks(:, :, sliderslice-ni+1));
end
set(handles.axes_outputimage, 'XTick', []);
set(handles.axes_outputimage, 'YTick', []);

% Show the slice number in text_outputslicenumber.
set(handles.text_outputslicenumber, 'String', sliderslice);

% --- Executes during object creation, after setting all properties.
function slider_outputslices_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% ________________________ INPUTS and TEXT BOXES __________________________
% _________________________________________________________________________
function text_imagename_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_imagename_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function text_inputslicenumber_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_inputslicenumber_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function text_outputslicenumber_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_outputslicenumber_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function listbox_log_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function listbox_log_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popup_plane_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function popup_plane_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function text_xres_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_xres_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function text_yres_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_yres_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function text_slicethickness_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_slicethickness_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function text_firstslice_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_firstslice_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function text_lastslice_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_lastslice_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function text_scalefactor_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_scalefactor_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function text_edgethreshold_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_edgethreshold_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function text_edgesigma_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_edgesigma_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function text_xdimension_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_xdimension_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function text_ydimension_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_ydimension_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function text_volume_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_volume_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function togglebutton_action_Callback(hObject, eventdata, handles)

function popup_show_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function popup_show_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

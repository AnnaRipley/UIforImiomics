function varargout = del1(varargin)
%DEL1 MATLAB code file for del1.fig
%      DEL1, by itself, creates a new DEL1 or raises the existing
%      singleton*.
%
%      H = DEL1 returns the handle to a new DEL1 or the handle to
%      the existing singleton*.
%
%      DEL1('Property','Value',...) creates a new DEL1 using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to del1_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      DEL1('CALLBACK') and DEL1('CALLBACK',hObject,...) call the
%      local function named CALLBACK in DEL1.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help del1

% Last Modified by GUIDE v2.5 31-May-2018 13:01:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @del1_OpeningFcn, ...
                   'gui_OutputFcn',  @del1_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before del1 is made visible.
function del1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Set to 1 when doing the user tests, 0 therwise
usertest = 0;

% Load pre-packaged data
if usertest
    %del1_script_usertest
else
    del1_script
end

% Volumes and other data are defined in the script, save handles when needed.
handles.corrVol = corrVol;
handles.corrdata = datamatrix;
handles.noCorrMaps = noFiles;

handles.imistruct = imiFatVolarray;
sizearray = size(imiFatVolarray);
handles.structsize = sizearray(2); % should be equal to poemlength

handles.female = 1; % true or false, depends on current Rvol
handles.dxavalsFemale = dxavalsFemale;
handles.dxavalsMale = dxavalsMale;
handles.noFemales = noFemales;
handles.noMales = noMales;
handles.IDsFemale = poemIDsFemale;
handles.IDsMale = poemIDsMale;

handles.xval = xinit;
handles.yval = yinit;
handles.zval = zinit;
handles.xinit = xinit;
handles.yinit = yinit;
handles.zinit = zinit;
handles.xmax = xmax;
handles.ymax = ymax;
handles.zmax = zmax;

handles.corscale = corscale;
handles.sagscale = sagscale;
handles.axiscale = axiscale;

handles.fatscale = fatscale;
handles.jdscale = jdscale;

% q = current index in the subject group (starts at female subject 1)
handles.qval = 1;
handles.qmax = noFemales; % changes to noMales when needed

% Zoom info
handles.zoomNr = 0;
handles.radius = radius;
handles.xrange = 1:xmax;
handles.yrange = 1:ymax;
handles.zrange = 1:zmax;

handles.allowDPCB = 1; % allow DragPoint CallBacks?
handles.sourceNr = 1; % 1: Correlation vols (8 different kinds!), 2: Image volume X
handles.datNr = 1; % 1: Fat, 2: Jacobian det

% Current "selected" view
% 1: Coronal, 2: Sagittal, 3: Axial
handles.viewNr = 1;

% Save data for del2
q = handles.qval;
poemid = poemIDsFemale(q, 1);
IDselected = poemid;
save('IDselected.mat', 'IDselected');

% Set filenames
for f = 1 : noFiles
    filenameString{f} = filenames(f);
end
set(handles.popupmenuFilename, 'String', filenameString);
set(handles.popupmenuFilename, 'Value', 4);
handles.currentCorr = 4;
handles.currentVol = squeeze( handles.corrVol(:,:,:,handles.currentCorr) );

% Set current Fat and Jacobian Det volumes
q = handles.qval;
% current Corrvol male / female ?
if handles.female
    ind = handles.IDsFemale(q, 2);
else
    ind = handles.IDsMale(q, 2);
end
handles.currentFat = double( handles.imistruct(ind).DefVol ) / handles.fatscale;
handles.currentJD = double( handles.imistruct(ind).JacDet ) / handles.jdscale;

% Custom colormap
handles.BWRmap = BWRmap;


% axesCor : coronal / front

handles.sparefCor = imref2d([ymax xmax]); % spatial reference
handles.sliceCor = squeeze( handles.currentVol(:,:,handles.zval) );
handles.hImageCor = imshow(handles.sliceCor', handles.sparefCor, 'Parent', handles.axesCor, 'Colormap', jet );
set(handles.axesCor, 'DataAspectRatio', corscale);
set(handles.hImageCor, 'ButtonDownFcn', @axesCor_ButtonDownFcn);
%handles.colbar = colorbar('southOutside');
handles.axesCor.XTick = [];
handles.axesCor.YTick = [];
handles.dragpointCor = new_dragpointCor(handles); % dragpoint


% axesSag : sagittal / side

handles.sparefSag = imref2d([handles.ymax handles.zmax]);
handles.sliceSag = squeeze( handles.currentVol(handles.xval,:,:) );
handles.hImageSag = imshow(handles.sliceSag, handles.sparefSag, 'Parent', handles.axesSag, 'Colormap', jet );
set(handles.axesSag, 'DataAspectRatio', sagscale);
set(handles.hImageSag, 'ButtonDownFcn', @axesSag_ButtonDownFcn);
handles.axesSag.XTick = [];
handles.axesSag.YTick = [];
handles.dragpointSag = new_dragpointSag(handles); % dragpoint
handles.colbar = colorbar(handles.axesSag, 'westOutside');


% axesAxi : axial / transverse

handles.sparefAxi = imref2d([handles.xmax handles.zmax]);
handles.sliceAxi = squeeze( handles.currentVol(:,handles.yval,:) );
handles.hImageAxi = imshow(handles.sliceAxi', handles.sparefAxi, 'Parent', handles.axesAxi, 'Colormap', jet );
set(handles.axesAxi, 'DataAspectRatio', axiscale);
set(handles.hImageAxi, 'ButtonDownFcn', @axesAxi_ButtonDownFcn);
handles.axesAxi.XTick = [];
handles.axesAxi.YTick = [];
handles.dragpointAxi = new_dragpointAxi(handles); % dragpoint


% For the colorbar and scaling
set(handles.axesCor, 'CLim', [-1 1]);
set(handles.axesSag, 'CLim', [-1 1]);
set(handles.axesAxi, 'CLim', [-1 1]);


% axesDat
handles.rowFemale = startrowFemale;
handles.rowMale = startrowMale;
handles.Fatlims = [0 1];
handles.JDlims = [0 2];


% textboxes
set(handles.editX, 'String', int2str(handles.xval));
set(handles.editY, 'String', int2str(handles.yval));
set(handles.editZ, 'String', int2str(handles.zval));


% Choose default command line output for del1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
update_axesDat();
% UIWAIT makes del1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = del1_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function axesCor_CreateFcn(~, eventdata, handles)
% hObject    handle to axesCor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axesCor


% --- Executes during object creation, after setting all properties.
function axesDat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axesDat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axesDat


% X X X


function editX_Callback(hObject, eventdata, handles)
% hObject    handle to editX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
new_xval = round(str2double(get(hObject, 'String')));

if new_xval > 0 && new_xval <= handles.xmax
    handles.xval = new_xval;
    
    guidata(hObject, handles);
    update_axesSag();
    update_dragpointCor();
    update_dragpointAxi();
    update_axesDat();
else
    set(hObject, 'String', int2str(handles.xval));
end
% Hints: get(hObject,'String') returns contents of editX as text
%        str2double(get(hObject,'String')) returns contents of editX as a double


% --- Executes during object creation, after setting all properties.
function editX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonXup.
function pushbuttonXup_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonXup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x = handles.xval + 1;
handles.xval = coordlimit(x, handles.xmax);
set(handles.editX, 'String', int2str(handles.xval));

guidata(hObject, handles);
update_axesSag();
update_dragpointCor();
update_dragpointAxi();
update_axesDat();



% --- Executes on button press in pushbuttonXdown.
function pushbuttonXdown_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonXdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x = handles.xval - 1;
handles.xval = coordlimit(x, handles.xmax); 
set(handles.editX, 'String', int2str(handles.xval));

guidata(hObject, handles);
update_axesSag();
update_dragpointCor();
update_dragpointAxi();
update_axesDat();



% Y Y Y


function editY_Callback(hObject, eventdata, handles)
% hObject    handle to editY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
new_yval = round(str2double(get(hObject, 'String')));

if new_yval > 0 && new_yval <= handles.ymax
    handles.yval = new_yval;

    guidata(hObject, handles);
    update_axesAxi();
    update_dragpointCor();
    update_dragpointSag();
    update_axesDat();
else
    set(hObject, 'String', int2str(handles.yval));
end
% Hints: get(hObject,'String') returns contents of editY as text
%        str2double(get(hObject,'String')) returns contents of editY as a double


% --- Executes during object creation, after setting all properties.
function editY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonYup.
function pushbuttonYup_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonYup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
y = handles.yval + 1;
handles.yval = coordlimit(y, handles.ymax);
set(handles.editY, 'String', int2str(handles.yval));

guidata(hObject, handles);
update_axesAxi();
update_dragpointCor();
update_dragpointSag();
update_axesDat();



% --- Executes on button press in pushbuttonYdown.
function pushbuttonYdown_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonYdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
y = handles.yval - 1;
handles.yval = coordlimit(y, handles.ymax);
set(handles.editY, 'String', int2str(handles.yval));

guidata(hObject, handles);
update_axesAxi();
update_dragpointCor();
update_dragpointSag();
update_axesDat();



% Z Z Z


function editZ_Callback(hObject, eventdata, handles)
% hObject    handle to editZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
new_zval = round(str2double(get(hObject, 'String')));

if new_zval > 0 && new_zval <= handles.zmax
    handles.zval = new_zval;

    guidata(hObject, handles);

    update_axesCor();
    update_dragpointSag();
    update_dragpointAxi();
    update_axesDat();
else    
    set(hObject, 'String', int2str(handles.zval));
end
% Hints: get(hObject,'String') returns contents of editZ as text
%        str2double(get(hObject,'String')) returns contents of editZ as a double


% --- Executes during object creation, after setting all properties.
function editZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonZup.
function pushbuttonZup_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonZup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
z = handles.zval + 1;
handles.zval = coordlimit(z, handles.zmax);
set(handles.editZ, 'String', int2str(handles.zval));

guidata(hObject, handles);
update_axesCor();
update_dragpointSag();
update_dragpointAxi();
update_axesDat();



% --- Executes on button press in pushbuttonZdown.
function pushbuttonZdown_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonZdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
z = handles.zval - 1;
handles.zval = coordlimit(z, handles.zmax);
set(handles.editZ, 'String', int2str(handles.zval));

guidata(hObject, handles);
update_axesCor();
update_dragpointSag();
update_dragpointAxi();
update_axesDat();



% Update Images

% Updates the image in axesCor - if zval changed
function update_axesCor()
handles = guidata(gcf);

% Needed if zoomed
xr = handles.xrange;
yr = handles.yrange;

switch handles.sourceNr
    case 1
    handles.sliceCor = squeeze( handles.currentVol(xr,yr,handles.zval) );
    case 2
    switch handles.datNr
        case 1
        handles.sliceCor = squeeze( handles.currentFat(xr,yr,handles.zval) );
        case 2
        handles.sliceCor = squeeze( handles.currentJD(xr,yr,handles.zval) );
    end
end

set(handles.hImageCor, 'CData', handles.sliceCor');

new_z = handles.zval; 
value = (new_z - 1) / handles.zmax; % z min = 1 
set(handles.sliderCor, 'Value', value);

guidata(gcf, handles);
update_dragpointCor();

% Updates the image in axesSag - if xval changed
function update_axesSag()
handles = guidata(gcf);

% Needed if zoomed
yr = handles.yrange;
zr = handles.zrange;

switch handles.sourceNr
    case 1
    handles.sliceSag = squeeze( handles.currentVol(handles.xval,yr,zr) );
    case 2
    switch handles.datNr
        case 1
        handles.sliceSag = squeeze( handles.currentFat(handles.xval,yr,zr) );
        case 2
        handles.sliceSag = squeeze( handles.currentJD(handles.xval,yr,zr) );
    end
end
set(handles.hImageSag, 'CData', handles.sliceSag);

new_x = handles.xval; 
value = (new_x - 1) / handles.xmax; % x min = 1 
set(handles.sliderSag, 'Value', value);

guidata(gcf, handles)
update_dragpointSag();

% Updates the image in axesAxi - if yval changed
function update_axesAxi()
handles = guidata(gcf);

% Needed if zoomed
xr = handles.xrange;
zr = handles.zrange;

switch handles.sourceNr
    case 1
    handles.sliceAxi = squeeze( handles.currentVol(xr,handles.yval,zr) );
    case 2
    switch handles.datNr
        case 1
        handles.sliceAxi = squeeze( handles.currentFat(xr,handles.yval,zr) );
        case 2
        handles.sliceAxi = squeeze( handles.currentJD(xr,handles.yval,zr) );
    end
end
set(handles.hImageAxi, 'CData', handles.sliceAxi');

new_y = handles.yval; 
value = (new_y - 1) / handles.ymax; % y min = 1 
set(handles.sliderAxi, 'Value', value);

guidata(gcf, handles);
update_dragpointAxi();

% Updates graph in axesDat
function update_axesDat()

%sprintf("Running update_axesDat()")
handles = guidata(gcf);

imax = handles.structsize;
q = handles.qval;
x = handles.xval;
y = handles.yval;
z = handles.zval;

fi = 1;
mi = 1;

% Pick out data from voxel (x,y,z) of all male/female volumes
switch handles.datNr
    case 1 % Fat
        %sprintf("datNr = 1 : Fat") 
        JD = 0;
        for i = 1:imax
            if handles.imistruct(i).female
                handles.rowFemale(fi) = handles.imistruct(i).DefVol(x,y,z);
                fi = fi + 1;
            else
                handles.rowMale(mi) = handles.imistruct(i).DefVol(x,y,z);
                mi = mi + 1;
            end
        end
    case 2 % Jacobian Det
        %sprintf("datNr = 2 : Jacobian")
        JD = 1;
        for i = 1:imax
            if handles.imistruct(i).female
                handles.rowFemale(fi) = handles.imistruct(i).JacDet(x,y,z);
                fi = fi + 1;
            else
                handles.rowMale(mi) = handles.imistruct(i).JacDet(x,y,z);
                mi = mi + 1;
            end
        end
end

% Prepare for plotting
if handles.female
    dxavals = handles.dxavalsFemale;
    maxdxa = max(dxavals);
    if JD
        rowdouble = double(handles.rowFemale) / handles.jdscale;       
    else
        rowdouble = double(handles.rowFemale) / handles.fatscale;
    end
    handles.rowdouble = rowdouble;
    rowmax = max(rowdouble);
    noFemales = handles.noFemales;

    baseX = zeros(noFemales, 1);
    baseY = zeros(noFemales, 1);    
    for i = 1:noFemales
        baseX(i) = i / noFemales * maxdxa;
        baseY(i) = i / noFemales * rowmax; 
    end   
else
    dxavals = handles.dxavalsMale;
    maxdxa = max(dxavals);
    if JD
        rowdouble = double(handles.rowMale) / handles.jdscale;       
    else
        rowdouble = double(handles.rowMale) / handles.fatscale;
    end
    handles.rowdouble = rowdouble;
    rowmax = max(rowdouble);
    noMales = handles.noMales;
        
    baseX = zeros(noMales, 1);
    baseY = zeros(noMales, 1);    
    for i = 1:noMales
        baseX(i) = i / noMales * maxdxa;
        baseY(i) = i / noMales * rowmax; 
    end
end

% Find the best possible fitted line
p = polyfit(dxavals, rowdouble, 1);
polyY = polyval(p, baseX);

handles.hDataPlot = plot(handles.axesDat, dxavals, rowdouble, 'k .', dxavals(q), rowdouble(q), 'b square' , baseX, polyY, 'g -');  

% To make it clickable (the plot, not the axis)
set(handles.hDataPlot, 'ButtonDownFcn', @axesDat_ButtonDownFcn);

title(handles.axesDat, 'Scatter plot of data points')
xlabel(handles.axesDat, 'DXA values');

% Set limits for Fat / Jacobian
switch handles.datNr
    case 1
        ylabel(handles.axesDat, 'Fat fraction (deformed)');
        ylim(handles.axesDat, handles.Fatlims);
    case 2
        ylabel(handles.axesDat, 'Jacobian determinant');
        ylim(handles.axesDat, handles.JDlims);
end

guidata(gcf, handles);

update_textValue(handles)


% Dragpoint Callbacks

% --- Exectes when dragpoint is moved (in axesCor)
% pos = [x y]
function dragpointCallbackCor(pos)
handles = guidata(gcf);
znr = handles.zoomNr;
xrange = handles.xrange;
yrange = handles.yrange;
if handles.allowDPCB
    %sprintf("Callback Cor running")
    if znr == 0
        x = round(pos(1));
        y = round(pos(2));
    else
        xlen = size(xrange);
        xplace = pos(1)/handles.xmax * xlen(2);
        xind = max(round(xplace), 1);
        x = xrange(xind);
        ylen = size(yrange);
        yplace = pos(2)/handles.ymax * ylen(2);
        yind = max(round(yplace), 1);
        y = yrange(yind);
    end             
    handles.xval = coordlimit(x, handles.xmax);
    handles.yval = coordlimit(y, handles.ymax);
    set(handles.editX, 'String', int2str(x));
    set(handles.editY, 'String', int2str(y));
    
    guidata(gcf, handles);    
    update_axesSag();
    update_axesAxi();
    update_axesDat();
else
    % no action
    %sprintf("Callback Cor disallowed")
end

% pos = [z y]
function dragpointCallbackSag(pos)
handles = guidata(gcf);
znr = handles.zoomNr;
yrange = handles.yrange;
zrange = handles.zrange;
if handles.allowDPCB
    %sprintf("Callback Sag running")
    if znr == 0
        z = round(pos(1));    
        y = round(pos(2));
    else
        zlen = size(zrange);
        zplace = pos(1)/handles.zmax * zlen(2);
        zind = max(round(zplace), 1);
        z = zrange(zind);
        ylen = size(yrange);
        yplace = pos(2)/handles.ymax * ylen(2);
        yind = max(round(yplace), 1);
        y = yrange(yind);
    end
    handles.zval = coordlimit(z, handles.zmax);
    handles.yval = coordlimit(y, handles.ymax);
    set(handles.editZ, 'String', int2str(z));
    set(handles.editY, 'String', int2str(y));
    
    guidata(gcf, handles);    
    update_axesCor();
    update_axesAxi();
    update_axesDat();
else
    % no action
    %sprintf("Callback Sag disallowed")
end

% pos = [x z]
function dragpointCallbackAxi(pos)
handles = guidata(gcf);
znr = handles.zoomNr;
xrange = handles.xrange;
zrange = handles.zrange;
if handles.allowDPCB
    %sprintf("Callback Axi running")
    if znr == 0
        x = round(pos(1));
        z = round(pos(2));
    else
        xlen = size(xrange);
        xplace = pos(1)/handles.xmax * xlen(2);
        xind = max(round(xplace), 1);
        x = xrange(xind);
        zlen = size(zrange);
        zplace = pos(2)/handles.zmax * zlen(2);
        zind = max(round(zplace), 1);
        z = zrange(zind);
    end
    handles.xval = coordlimit(x, handles.xmax);
    handles.zval = coordlimit(z, handles.zmax); 
    set(handles.editX, 'String', int2str(x));
    set(handles.editZ, 'String', int2str(z));
    
    guidata(gcf, handles);    
    update_axesCor();
    update_axesSag();
    update_axesDat();
else
    % no action
    %sprintf("Callback Axi disallowed")
end


% New and Update Dragpoints


function dp = new_dragpointCor(handles)
dp = impoint(handles.axesCor, handles.xval, handles.yval);
setColor(dp, 'black'); %black
dpCor = @dragpointCallbackCor;
addNewPositionCallback(dp, dpCor);

function update_dragpointCor()
handles = guidata(gcf);
handles.allowDPCB = 0;
guidata(gcf, handles);
znr = handles.zoomNr;
xval = handles.xval;
yval = handles.yval;
if znr == 0
    % keep these
else
    [isit, ind] = ismember(xval, handles.xrange);
    len = size(handles.xrange);
    kvot = ind / len(2); 
    xval = kvot * handles.xmax;
    [isit, ind] = ismember(yval, handles.yrange);
    len = size(handles.yrange);
    kvot = ind / len(2); 
    yval = kvot * handles.ymax;
end

setPosition(handles.dragpointCor, xval, yval);
handles.allowDPCB = 1;
guidata(gcf, handles);

function dp = new_dragpointSag(handles)
dp = impoint(handles.axesSag, handles.zval, handles.yval);
setColor(dp, 'black');
dpSag = @dragpointCallbackSag;
addNewPositionCallback(dp, dpSag);

function update_dragpointSag()
handles = guidata(gcf);
handles.allowDPCB = 0;
guidata(gcf, handles);
znr = handles.zoomNr;
zval = handles.zval;
yval = handles.yval;
if znr == 0
    % keep these
else
    [isit, ind] = ismember(zval, handles.zrange);
    len = size(handles.zrange);
    kvot = ind / len(2); 
    zval = kvot * handles.zmax;
    [isit, ind] = ismember(yval, handles.yrange);
    len = size(handles.yrange);
    kvot = ind / len(2); 
    yval = kvot * handles.ymax;
end
setPosition(handles.dragpointSag, zval, yval);
handles.allowDPCB = 1;
guidata(gcf, handles);

function dp = new_dragpointAxi(handles) 
dp = impoint(handles.axesAxi, handles.xval, handles.zval);
setColor(dp, 'black');
dpAxi = @dragpointCallbackAxi;
addNewPositionCallback(dp, dpAxi);

function update_dragpointAxi()
handles = guidata(gcf);
handles.allowDPCB = 0;
guidata(gcf, handles);
znr = handles.zoomNr;
zval = handles.zval;
xval = handles.xval;
if znr == 0
    % keep these
else
    [isit, ind] = ismember(zval, handles.zrange);
    len = size(handles.zrange);
    kvot = ind / len(2); 
    zval = kvot * handles.zmax;
    [isit, ind] = ismember(xval, handles.xrange);
    len = size(handles.xrange);
    kvot = ind / len(2); 
    xval = kvot * handles.xmax;
end
setPosition(handles.dragpointAxi, xval, zval);
handles.allowDPCB = 1;
guidata(gcf, handles);



% Q Q Q


function editQ_Callback(hObject, eventdata, handles)
% hObject    handle to editQ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
new_qval = round(str2double(get(hObject, 'String')));
if new_qval > 0 && new_qval <= handles.qmax
    handles.qval = new_qval;
    guidata(hObject, handles);
    
    after_update_q(handles);
else
    % no change
    set(hObject, 'String', int2str(handles.qval));
end
% Hints: get(hObject,'String') returns contents of editQ as text
%        str2double(get(hObject,'String')) returns contents of editQ as a double


% --- Executes during object creation, after setting all properties.
function editQ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editQ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonQup.
function pushbuttonQup_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonQup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.qval < handles.qmax
    new_qval = handles.qval + 1;
    handles.qval = new_qval;
    set(handles.editQ, 'String', int2str(handles.qval));
    guidata(hObject, handles);
    
    after_update_q(handles);
else
    %no change
    handles.qval = handles.qmax;
    set(handles.editQ, 'String', int2str(handles.qval));
    guidata(hObject, handles);
end


% --- Executes on button press in pushbuttonQdown.
function pushbuttonQdown_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonQdown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.qval > 1
    new_qval = handles.qval - 1;
    handles.qval = new_qval;
    set(handles.editQ, 'String', int2str(handles.qval));
    guidata(hObject, handles);
    
    after_update_q(handles);
else
    %no change
    handles.qval = 1;
    set(handles.editQ, 'String', int2str(handles.qval));
    guidata(hObject, handles);
end


function after_update_q(handles)
q = handles.qval;

if handles.female
    poemid = handles.IDsFemale(q, 1);
    ind = handles.IDsFemale(q, 2);
else
    poemid = handles.IDsMale(q, 1);
    ind = handles.IDsMale(q, 2);
end

% Save poemid
set(handles.textIDnr, 'String', int2str(poemid) );
IDselected = poemid;
save('IDselected.mat', 'IDselected');
    
% Set new current image vols
handles.currentFat = double( handles.imistruct(ind).DefVol ) / handles.fatscale;
handles.currentJD = double( handles.imistruct(ind).JacDet ) / handles.jdscale;

guidata(gcf, handles);
    
update_axesDat();
if handles.sourceNr == 2
    update_axesCor();
    update_axesSag();
    update_axesAxi();
end
    

% Popup menus


% --- Executes on selection change in popupmenuSource.
function popupmenuSource_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuSource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(hObject, 'Value');
handles.sourceNr = value;

corrVal = handles.currentCorr;
pval = handles.corrdata(corrVal, 2);

if handles.corrdata(corrVal, 3)
    handles.datNr = 1; % 1: Fat, 2: Jacobian det
else
    handles.datNr = 2;
end

guidata(hObject, handles);

switch value
    case 1
        % Correlation Map
        colormap(handles.axesCor, jet );
        colormap(handles.axesSag, jet );
        colormap(handles.axesAxi, jet );
        set(handles.popupmenuColormap, 'Value', 1);
        if pval
            set(handles.axesCor, 'CLim', [0 1]);
            set(handles.axesSag, 'CLim', [0 1]);
            set(handles.axesAxi, 'CLim', [0 1]);
        else
            set(handles.axesCor, 'CLim', [-1 1]);
            set(handles.axesSag, 'CLim', [-1 1]);
            set(handles.axesAxi, 'CLim', [-1 1]);
        end

        setColor(handles.dragpointCor, 'black');
        setColor(handles.dragpointSag, 'black');
        setColor(handles.dragpointAxi, 'black');
        set(handles.popupmenuColor, 'Value', 1);
        
        guidata(gcf, handles);
        
        update_axesCor();
        update_axesSag();
        update_axesAxi();
        
    case 2
        % Specific Image Volume

        
        if handles.datNr == 1 % fat
            Fatlims = handles.Fatlims; % should be [0 1]
            set(handles.axesCor, 'CLim', Fatlims);
            set(handles.axesSag, 'CLim', Fatlims);
            set(handles.axesAxi, 'CLim', Fatlims);
            colormap(handles.axesCor, gray );
            colormap(handles.axesSag, gray );
            colormap(handles.axesAxi, gray );
            set(handles.popupmenuColormap, 'Value', 2);
            setColor(handles.dragpointCor, 'magenta');
            setColor(handles.dragpointSag, 'magenta');
            setColor(handles.dragpointAxi, 'magenta');
            set(handles.popupmenuColor, 'Value', 3);
        else
            % datNr = 2 = Jacobian det [0 2]
            JDlims = handles.JDlims;
            set(handles.axesCor, 'CLim', JDlims);
            set(handles.axesSag, 'CLim', JDlims);
            set(handles.axesAxi, 'CLim', JDlims);
            colormap(handles.axesCor, handles.BWRmap );
            colormap(handles.axesSag, handles.BWRmap );
            colormap(handles.axesAxi, handles.BWRmap );
            set(handles.popupmenuColormap, 'Value', 5);
            setColor(handles.dragpointCor, 'black');
            setColor(handles.dragpointSag, 'black');
            setColor(handles.dragpointAxi, 'black');
            set(handles.popupmenuColor, 'Value', 1);
        end

        guidata(gcf, handles);
        
        update_axesCor();
        update_axesSag();
        update_axesAxi();
        
    otherwise
        sprintf("Invalid value of popupmenuSource!");

end
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuSource contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuSource


% --- Executes during object creation, after setting all properties.
function popupmenuSource_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuSource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuFilename.
function popupmenuFilename_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuFilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = get(hObject, 'Value');
handles.currentCorr = val;
handles.currentVol = squeeze( handles.corrVol(:,:,:,val) );

guidata(gcf, handles);

fem = handles.corrdata(val, 1);
handles.female = fem;
if fem
    radiobuttonFemale_Callback(handles.radiobuttonFemale, eventdata, handles);
    %sprintf("female")
else
    radiobuttonMale_Callback(handles.radiobuttonMale, eventdata, handles);
    %sprintf("male")
end

%handles.pval ?
pval = handles.corrdata(val, 2);
if pval
    radiobuttonP_Callback(handles.radiobuttonP, eventdata, handles);
else
    radiobuttonR_Callback(handles.radiobuttonR, eventdata, handles);
end

datval = handles.corrdata(val, 3);
if datval % fat = 1
    radiobuttonFat_Callback(handles.radiobuttonFat, eventdata, handles);
else % JacDet
    radiobuttonJD_Callback(handles.radiobuttonJD, eventdata, handles);
end
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuFilename contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuFilename


% --- Executes during object creation, after setting all properties.
function popupmenuFilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuFilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuColor.
function popupmenuColor_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuColor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Marker color
% 1: Black, 2: White 3: Magenta
value = get(hObject, 'Value'); 
switch value
    case 1
    setColor(handles.dragpointCor, 'black');
    setColor(handles.dragpointSag, 'black');
    setColor(handles.dragpointAxi, 'black');    
    case 2
    setColor(handles.dragpointCor, 'white');
    setColor(handles.dragpointSag, 'white');
    setColor(handles.dragpointAxi, 'white');
    case 3
    setColor(handles.dragpointCor, 'magenta');
    setColor(handles.dragpointSag, 'magenta');
    setColor(handles.dragpointAxi, 'magenta');
end
guidata(gcf, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuColor contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuColor


% --- Executes during object creation, after setting all properties.
function popupmenuColor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuColor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuColormap.
function popupmenuColormap_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuColormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(hObject,'Value');
switch value
    case 1
        colormap(handles.axesCor, jet );
        colormap(handles.axesSag, jet );
        colormap(handles.axesAxi, jet );
    case 2
        colormap(handles.axesCor, gray );
        colormap(handles.axesSag, gray );
        colormap(handles.axesAxi, gray );
    case 3
        colormap(handles.axesCor, parula );
        colormap(handles.axesSag, parula );
        colormap(handles.axesAxi, parula );
    case 4
        colormap(handles.axesCor, hot );
        colormap(handles.axesSag, hot );
        colormap(handles.axesAxi, hot );
    case 5
        colormap(handles.axesCor, handles.BWRmap );
        colormap(handles.axesSag, handles.BWRmap );
        colormap(handles.axesAxi, handles.BWRmap );
end
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuColormap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuColormap


% --- Executes during object creation, after setting all properties.
function popupmenuColormap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuColormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% Key shortcuts

% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
keyPressed = eventdata.Key;
switch keyPressed
    case 'leftarrow'
        uicontrol(handles.pushbuttonQdown);
        pushbuttonQdown_Callback(hObject, eventdata, handles);
    case 'rightarrow'
        uicontrol(handles.pushbuttonQup);
        pushbuttonQup_Callback(hObject, eventdata, handles);
    case 'x'
        uicontrol(handles.editX);
        handles.viewNr = 2; % sag
        guidata(hObject, handles);
    case 'y'
        uicontrol(handles.editY);
        handles.viewNr = 3; % axi
        guidata(hObject, handles);
    case 'z'
        uicontrol(handles.editZ);
        handles.viewNr = 1; % cor
        guidata(hObject, handles);
    case 'r'
        uicontrol(handles.pushbuttonReset);
        pushbuttonReset_Callback(hObject, eventdata, handles);
    case 'i'
        uicontrol(handles.pushbuttonZoomIn);
        pushbuttonZoomIn_Callback(hObject, eventdata, handles);
    case 'o'
        uicontrol(handles.pushbuttonZoomOut);
        pushbuttonZoomOut_Callback(hObject, eventdata, handles);
    case 'f'
        uicontrol(handles.pushbuttonAutoLims);
        pushbuttonAutoLims_Callback(hObject, eventdata, handles);
        
    otherwise
        %sprintf("That key is unbound.")
end


% Select Axes

% --- Executes on mouse press over axes background.
function axesCor_ButtonDownFcn(hObject, ~, ~)
% hObject    handle to axesCor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
handles.viewNr = 1;
set(handles.textScrollCor, 'Visible', 'on');
set(handles.textScrollSag, 'Visible', 'off');
set(handles.textScrollAxi, 'Visible', 'off');
%sprintf("Cor button down")

point = get(handles.axesCor, 'CurrentPoint');
imX = point(1, 1);
imY = point(1, 2);

guidata(gcf, handles);

setPosition(handles.dragpointCor, imX, imY);

uicontrol(handles.editZ);


% --- Executes on mouse press over axes background.
function axesSag_ButtonDownFcn(hObject, ~, ~)
% hObject    handle to axesSag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
handles.viewNr = 2;
set(handles.textScrollCor, 'Visible', 'off');
set(handles.textScrollSag, 'Visible', 'on');
set(handles.textScrollAxi, 'Visible', 'off');
%sprintf("Sag button down")

point = get(handles.axesSag, 'CurrentPoint');
imX = point(1, 1);
imY = point(1, 2);

guidata(gcf, handles);

setPosition(handles.dragpointSag, imX, imY);

uicontrol(handles.editX);


% --- Executes on mouse press over axes background.
function axesAxi_ButtonDownFcn(hObject, ~, ~)
% hObject    handle to axesAxi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
handles.viewNr = 3;
set(handles.textScrollCor, 'Visible', 'off');
set(handles.textScrollSag, 'Visible', 'off');
set(handles.textScrollAxi, 'Visible', 'on');
%sprintf("Axi button down")

point = get(handles.axesAxi, 'CurrentPoint');
imX = point(1, 1);
imY = point(1, 2);

guidata(gcf, handles);

setPosition(handles.dragpointAxi, imX, imY);

uicontrol(handles.editY);


% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)
scroll = eventdata.VerticalScrollCount;

if scroll < 0
    switch handles.viewNr
        case 1
            uicontrol(handles.pushbuttonZup);
            pushbuttonZup_Callback(hObject, eventdata, handles);
        case 2
            uicontrol(handles.pushbuttonXup);
            pushbuttonXup_Callback(hObject, eventdata, handles);
        case 3
            uicontrol(handles.pushbuttonYup);
            pushbuttonYup_Callback(hObject, eventdata, handles);
    end
else
    switch handles.viewNr
        case 1
            uicontrol(handles.pushbuttonZdown);
            pushbuttonZdown_Callback(hObject, eventdata, handles);
        case 2
            uicontrol(handles.pushbuttonXdown);
            pushbuttonXdown_Callback(hObject, eventdata, handles);
        case 3
            uicontrol(handles.pushbuttonYdown);
            pushbuttonYdown_Callback(hObject, eventdata, handles);
    end
end


function update_textValue(handles)
x = handles.xval;
y = handles.yval;
z = handles.zval;
q = handles.qval;

imvol_value = handles.rowdouble(q);

set(handles.textValue, 'String', num2str(imvol_value) );

corrNr = handles.currentCorr;
if handles.corrdata(corrNr, 2) % pval
    Rval = handles.corrVol(x,y,z,corrNr + 2);
    Pval = handles.corrVol(x,y,z,corrNr);
else
    Rval = handles.corrVol(x,y,z,corrNr);
    Pval = handles.corrVol(x,y,z,corrNr - 2);
end

set(handles.textCorVal, 'String', num2str(Rval) );
set(handles.textPval, 'String', num2str(Pval) );

guidata(gcf, handles);


% Sliders


% --- Executes on slider movement.
function sliderSag_Callback(hObject, eventdata, handles)
% hObject    handle to sliderSag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(hObject, 'Value');

new_x = max(round(handles.xmax * value), 1); % x min = 1
handles.xval = new_x;

set(handles.editX, 'String', int2str(handles.xval));
guidata(hObject, handles);
% update axes
update_axesSag();
update_dragpointCor();
update_dragpointAxi();
update_axesDat();
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sliderSag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderSag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderCor_Callback(hObject, eventdata, handles)
% hObject    handle to sliderCor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(hObject, 'Value');

new_z = max(round(handles.zmax * value), 1); % z min = 1
handles.zval = new_z;

set(handles.editZ, 'String', int2str(handles.zval));
guidata(hObject, handles);
% update axes
update_axesCor();
update_dragpointSag();
update_dragpointAxi();
update_axesDat();
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sliderCor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderCor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderAxi_Callback(hObject, eventdata, handles)
% hObject    handle to sliderAxi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(hObject, 'Value');

new_y = max(round(handles.ymax * value), 1); % y min = 1
handles.yval = new_y;

set(handles.editY, 'String', int2str(handles.yval));
guidata(hObject, handles);
% update axes
update_axesAxi();
update_dragpointCor();
update_dragpointSag();
update_axesDat();
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sliderAxi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderAxi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbuttonReset.
function pushbuttonReset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.xval = handles.xinit;
handles.yval = handles.yinit;
handles.zval = handles.zinit;
set(handles.editX, 'String', int2str(handles.xval));
set(handles.editY, 'String', int2str(handles.yval));
set(handles.editZ, 'String', int2str(handles.zval));
% Reset Zoom
handles.xrange = 1:handles.xmax;
handles.yrange = 1:handles.ymax;
handles.zrange = 1:handles.zmax;
handles.zoomNr = 0;
guidata(hObject, handles);
update_axesCor();
update_axesSag();
update_axesAxi();
update_axesDat();


% Limits x,y,z coords from being < 1 or > max 
function aout = coordlimit(a, alim)
aout = a;
if a < 1
    aout = 1;
end
if a > alim
    aout = alim;
end


% --- Executes on mouse press over axes background.
function axesDat_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesDat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject); % because need to be updated
point = get(handles.axesDat, 'CurrentPoint');
%sprintf("axes Data button down")
xdxa = point(1, 1);
fem = handles.female;
if fem
    dxavals = handles.dxavalsFemale;
else
    dxavals = handles.dxavalsMale;
end
differenceX = (dxavals - xdxa).^2;

Yscaleup = 10^5; % weight on diffY
yrow = Yscaleup * point(1, 2);
differenceY = (Yscaleup * handles.rowdouble - yrow).^2;
differenceTot = differenceX + differenceY; 
[val, ind] = min(differenceTot);

new_q = ind;

handles.qval = new_q; % new_q is within q limits
set(handles.editQ, 'String', int2str(new_q));
guidata(gcf, handles);
editQ_Callback(handles.editQ, eventdata, handles); % handles all other changes


% --- Executes on button press in pushbuttonZoomIn.
function pushbuttonZoomIn_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonZoomIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
maxzoom = 3;
r = handles.radius;
x = handles.xval;
y = handles.yval;
z = handles.zval;

znr = handles.zoomNr + 1;
if znr > maxzoom
    %do nothing
else
    handles.zoomNr = znr;
    %zoom in
    xlo = coordlimit(x - r *(4-znr), handles.xmax);
    xhi = coordlimit(x + r *(4-znr), handles.xmax);
    ylo = coordlimit(y - r *(4-znr), handles.ymax);
    yhi = coordlimit(y + r *(4-znr), handles.ymax);
    zlo = coordlimit(z - r *(4-znr), handles.zmax);
    zhi = coordlimit(z + r *(4-znr), handles.zmax);
    handles.xrange = xlo:xhi;
    handles.yrange = ylo:yhi;
    handles.zrange = zlo:zhi;
    
    guidata(hObject, handles);
    update_axesCor();
    update_axesSag();
    update_axesAxi();
end

% --- Executes on button press in pushbuttonZoomOut.
function pushbuttonZoomOut_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonZoomOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
znr = handles.zoomNr - 1;
if znr < 0
    %do nothing
else
    r = handles.radius;
    x = handles.xval;
    y = handles.yval;
    z = handles.zval;
    %zoom out
    if znr == 0
        handles.xrange = 1:handles.xmax;
        handles.yrange = 1:handles.ymax;
        handles.zrange = 1:handles.zmax;
    else
        xlo = coordlimit(x - r *(4-znr), handles.xmax);
        xhi = coordlimit(x + r *(4-znr), handles.xmax);
        ylo = coordlimit(y - r *(4-znr), handles.ymax);
        yhi = coordlimit(y + r *(4-znr), handles.ymax);
        zlo = coordlimit(z - r *(4-znr), handles.zmax);
        zhi = coordlimit(z + r *(4-znr), handles.zmax);
        handles.xrange = xlo:xhi;
        handles.yrange = ylo:yhi;
        handles.zrange = zlo:zhi;
    end
    handles.zoomNr = znr;
    
    guidata(hObject, handles);
    update_axesCor();
    update_axesSag();
    update_axesAxi();
end



% Data Ylim / Display range limits

function editLimUpper_Callback(hObject, eventdata, handles)
% hObject    handle to editLimUpper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = str2double(get(hObject, 'String'));
lowerupper = handles.JDlims;

if value > 5
    value = 5;
    set(hObject, 'String', int2str(value));
end
if value < 1
    value = 1;
    set(hObject, 'String', int2str(value));
end

if handles.sourceNr == 2
    set(handles.axesCor, 'CLim', [lowerupper(1) value]);
    set(handles.axesSag, 'CLim', [lowerupper(1) value]);
    set(handles.axesAxi, 'CLim', [lowerupper(1) value]);
end

set(handles.sliderLimUpper, 'Value', (value - 1) / 4);

handles.JDlims = [lowerupper(1) value];

guidata(hObject, handles);

if handles.sourceNr == 2
    update_axesCor();
    update_axesSag();
    update_axesAxi();
end
update_axesDat();
% Hints: get(hObject,'String') returns contents of editLimUpper as text
%        str2double(get(hObject,'String')) returns contents of editLimUpper as a double


% --- Executes during object creation, after setting all properties.
function editLimUpper_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLimUpper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderLimUpper_Callback(hObject, eventdata, handles)
% hObject    handle to sliderLimUpper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lowerupper = handles.JDlims;
% value [0 1]
value = get(hObject, 'Value');
new_upper = value * 4 + 1;

if handles.sourceNr == 2
    set(handles.axesCor, 'CLim', [lowerupper(1) new_upper]);
    set(handles.axesSag, 'CLim', [lowerupper(1) new_upper]);
    set(handles.axesAxi, 'CLim', [lowerupper(1) new_upper]);
end

set(handles.editLimUpper, 'String', num2str(new_upper));

handles.JDlims = [lowerupper(1) new_upper];

guidata(hObject, handles);

if handles.sourceNr == 2
    update_axesCor();
    update_axesSag();
    update_axesAxi();
end
update_axesDat();
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sliderLimUpper_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderLimUpper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function editLimLower_Callback(hObject, eventdata, handles)
% hObject    handle to editLimLower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = str2double(get(hObject, 'String'));
lowerupper = handles.JDlims;

if value < -3
    value = -3;
    set(hObject, 'String', int2str(value));
end
if value >= 1
    value = 0.99;
    set(hObject, 'String', int2str(value));
end

if handles.sourceNr == 2
    set(handles.axesCor, 'CLim', [value lowerupper(2)]);
    set(handles.axesSag, 'CLim', [value lowerupper(2)]);
    set(handles.axesAxi, 'CLim', [value lowerupper(2)]);
end

set(handles.sliderLimLower, 'Value', (value + 3) / 4);

handles.JDlims = [lowerupper(1) value];

guidata(hObject, handles);

if handles.sourceNr == 2
    update_axesCor();
    update_axesSag();
    update_axesAxi();
end
update_axesDat();
% Hints: get(hObject,'String') returns contents of editLimLower as text
%        str2double(get(hObject,'String')) returns contents of editLimLower as a double


% --- Executes during object creation, after setting all properties.
function editLimLower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLimLower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderLimLower_Callback(hObject, eventdata, handles)
% hObject    handle to sliderLimLower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lowerupper = handles.JDlims;
% value [0 1]
value = get(hObject, 'Value');
new_lower = (value - 1) * 4 + 1;

if new_lower > 0.99
    new_lower = 0.99;
end

if handles.sourceNr == 2
    set(handles.axesCor, 'CLim', [new_lower lowerupper(2)]);
    set(handles.axesSag, 'CLim', [new_lower lowerupper(2)]);
    set(handles.axesAxi, 'CLim', [new_lower lowerupper(2)]);
end

set(handles.editLimLower, 'String', num2str(new_lower));

handles.JDlims = [new_lower lowerupper(2)];

guidata(hObject, handles);

if handles.sourceNr == 2
    update_axesCor();
    update_axesSag();
    update_axesAxi();
end
update_axesDat();
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sliderLimLower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderLimLower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function update_limits(lower, upper)
handles = guidata(gcf);

if lower < -3
    lower = -3;
end
if lower > 0.99
    lower = 0.99;
end
if upper > 5
    upper = 5;
end
if upper < 1
    upper = 1;
end

set(handles.editLimLower, 'String', num2str(lower));
set(handles.editLimUpper, 'String', num2str(upper));
set(handles.sliderLimLower, 'Value', (lower + 3) / 4);
set(handles.sliderLimUpper, 'Value', (upper - 1) / 4);

handles.JDlims = [lower upper];

guidata(gcf, handles);

if handles.sourceNr == 2 && handles.datNr == 2
    set(handles.axesCor, 'CLim', [lower upper]);
    set(handles.axesSag, 'CLim', [lower upper]);
    set(handles.axesAxi, 'CLim', [lower upper]);
    
    update_axesCor();
    update_axesSag();
    update_axesAxi();
end
update_axesDat();


% --- Executes on button press in pushbuttonAutoLims.
function pushbuttonAutoLims_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAutoLims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
maxval = max(handles.rowdouble);
minval = min(handles.rowdouble);

distance = maxval - minval;

update_limits(minval - 0.05*distance, maxval + 0.05*distance);


% --- Executes on button press in pushbuttonResetLims.
function pushbuttonResetLims_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonResetLims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_limits(0, 2);


% Radio Buttons

% --- Executes on button press in radiobuttonJD.
function radiobuttonJD_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonJD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%sprintf("JD callback running")
handles.datNr = 2;
guidata(hObject, handles);

set(handles.radiobuttonJD, 'Value', 1);

oldNr = get(handles.popupmenuFilename, 'Value');
newrow = handles.corrdata(oldNr,:);
newrow(3) = 0; % JD
% update popupmenuFilename
newNr = set_new_filename(handles, newrow);
% update currentCorr
handles.currentCorr = newNr;

% update currentVol
handles.currentVol = squeeze( handles.corrVol(:,:,:,newNr) );

guidata(hObject, handles);

if handles.sourceNr == 2
    JDlims = handles.JDlims;
    set(handles.axesCor, 'CLim', JDlims);
    set(handles.axesSag, 'CLim', JDlims);
    set(handles.axesAxi, 'CLim', JDlims);
    
    colormap(handles.axesCor, handles.BWRmap );
    colormap(handles.axesSag, handles.BWRmap );
    colormap(handles.axesAxi, handles.BWRmap );
    set(handles.popupmenuColormap, 'Value', 5);
end

update_axesCor();
update_axesSag();
update_axesAxi();
update_axesDat();
% Hint: get(hObject,'Value') returns toggle state of radiobuttonJD


% --- Executes on button press in radiobuttonFat.
function radiobuttonFat_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonFat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%sprintf("Fat callback running")
handles.datNr = 1;
guidata(hObject, handles);

set(handles.radiobuttonFat, 'Value', 1);

oldNr = get(handles.popupmenuFilename, 'Value');
newrow = handles.corrdata(oldNr,:);
newrow(3) = 1; % Fat
% update popupmenuFilename
newNr = set_new_filename(handles, newrow);
% update currentCorr
handles.currentCorr = newNr;

% update currentVol
handles.currentVol = squeeze( handles.corrVol(:,:,:,newNr) );

guidata(hObject, handles);

if handles.sourceNr == 2
    Fatlims = handles.Fatlims;
    set(handles.axesCor, 'CLim', Fatlims);
    set(handles.axesSag, 'CLim', Fatlims);
    set(handles.axesAxi, 'CLim', Fatlims);
    
    colormap(handles.axesCor, gray );
    colormap(handles.axesSag, gray );
    colormap(handles.axesAxi, gray );
    set(handles.popupmenuColormap, 'Value', 2);
end

update_axesCor();
update_axesSag();
update_axesAxi();
update_axesDat();
% Hint: get(hObject,'Value') returns toggle state of radiobuttonFat


% --- Executes on button press in radiobuttonFemale.
function radiobuttonFemale_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonFemale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.radiobuttonFemale, 'Value', 1);
handles.female = 1;
handles.qmax = handles.noFemales;
if handles.qval > handles.qmax
    handles.qval = handles.qmax;
end
    
oldNr = get(handles.popupmenuFilename, 'Value');
newrow = handles.corrdata(oldNr,:);
newrow(1) = 1; % Female
% update popupmenuFilename
newNr = set_new_filename(handles, newrow);
% update currentCorr
handles.currentCorr = newNr;

% update currentVol
handles.currentVol = squeeze( handles.corrVol(:,:,:,newNr) );

guidata(hObject, handles);

after_update_q(handles);
% axesDat

update_axesCor();
update_axesSag();
update_axesAxi();
% Hint: get(hObject,'Value') returns toggle state of radiobuttonFemale


% --- Executes on button press in radiobuttonMale.
function radiobuttonMale_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonMale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.radiobuttonMale, 'Value', 1);
handles.female = 0;
handles.qmax = handles.noMales;
if handles.qval > handles.qmax
    handles.qval = handles.qmax;
end

oldNr = get(handles.popupmenuFilename, 'Value');
newrow = handles.corrdata(oldNr,:);
newrow(1) = 0; % Male
% update popupmenuFilename
newNr = set_new_filename(handles, newrow);
% update currentCorr
handles.currentCorr = newNr;

% update currentVol
handles.currentVol = squeeze( handles.corrVol(:,:,:,newNr) );

guidata(hObject, handles);

after_update_q(handles);
% axesDat

update_axesCor();
update_axesSag();
update_axesAxi();
% Hint: get(hObject,'Value') returns toggle state of radiobuttonMale


% --- Executes on button press in radiobuttonP.
function radiobuttonP_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.radiobuttonP, 'Value', 1);

oldNr = get(handles.popupmenuFilename, 'Value');
newrow = handles.corrdata(oldNr,:);
newrow(2) = 1; % P

if handles.sourceNr == 1
    set(handles.axesCor, 'CLim', [0 1]);
    set(handles.axesSag, 'CLim', [0 1]);
    set(handles.axesAxi, 'CLim', [0 1]);
end
% update popupmenuFilename
newNr = set_new_filename(handles, newrow);
% update currentCorr
handles.currentCorr = newNr;

% update currentVol
handles.currentVol = squeeze( handles.corrVol(:,:,:,newNr) );

guidata(hObject, handles);

update_axesCor();
update_axesSag();
update_axesAxi();
% Hint: get(hObject,'Value') returns toggle state of radiobuttonP


% --- Executes on button press in radiobuttonR.
function radiobuttonR_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.radiobuttonR, 'Value', 1);

oldNr = get(handles.popupmenuFilename, 'Value');
newrow = handles.corrdata(oldNr,:);
newrow(2) = 0; % R

if handles.sourceNr == 1
    set(handles.axesCor, 'CLim', [-1 1]);
    set(handles.axesSag, 'CLim', [-1 1]);
    set(handles.axesAxi, 'CLim', [-1 1]);
end
% update popupmenuFilename
newNr = set_new_filename(handles, newrow);
% update currentCorr
handles.currentCorr = newNr;

% update currentVol
handles.currentVol = squeeze( handles.corrVol(:,:,:,newNr) );

guidata(hObject, handles);

update_axesCor();
update_axesSag();
update_axesAxi();
% Hint: get(hObject,'Value') returns toggle state of radiobuttonR


function rowNr = set_new_filename(handles, newrow)
data = handles.corrdata;
maxRow = handles.noCorrMaps;
rowNr = 0;

for i = 1:maxRow
    if data(i,:) == newrow
        rowNr = i;
    end
end
if rowNr == 0
    sprintf("Error in set_new_filename()")
else
set(handles.popupmenuFilename, 'Value', rowNr);
end


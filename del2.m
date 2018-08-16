function varargout = del2(varargin)
% DEL2 MATLAB code for del2.fig
%      DEL2, by itself, creates a new DEL2 or raises the existing
%      singleton*.
%
%      H = DEL2 returns the handle to a new DEL2 or the handle to
%      the existing singleton*.
%
%      DEL2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DEL2.M with the given input arguments.
%
%      DEL2('Property','Value',...) creates a new DEL2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before del2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to del2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help del2

% Last Modified by GUIDE v2.5 31-May-2018 13:03:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @del2_OpeningFcn, ...
                   'gui_OutputFcn',  @del2_OutputFcn, ...
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


% --- Executes just before del2 is made visible.
function del2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to del2 (see VARARGIN)

% Set 1 for test, 0 otherwise
test = 0; % set to 1 to use a small test.mat instead of the .vtk files
handles.usertest = 0; % set to 1 if using the uidata_usertest.mat in part 1

handles.test = test;

if test
    %load('test.mat');
else
    del2_script
end

handles.pathDeformed = pathDeformed;
handles.pathUndeformed = pathUndeformed;

handles.poemID = IDselected;
handles.female = 1;
handles.poemdxa = POEMdxa;
handles.poemlength = poemlength;

handles.xmax = xmax;
handles.ymax = ymax;
handles.zmax = zmax;
handles.xval = xinit;
handles.yval = yinit;
handles.zval = zinit; 
handles.xinit = xinit;
handles.yinit = yinit;
handles.zinit = zinit; 
handles.corscale = corscale;
handles.sagscale = sagscale;
handles.axiscale = axiscale;

if test
    handles.jacdetvol = jacdetvol;
    handles.undefvol = undefvol;
    handles.deformvol = deformvol;
else
    handles.jacdetvol = zeros(xmax, ymax, zmax);
    handles.undefvol = zeros(xmax, ymax, zmax);
    handles.deformvol = zeros(xmax, ymax, zmax);
end

handles.refvolMale = refvolMale;
handles.refvolFemale = refvolFemale;

handles.diffvol = handles.deformvol - refvolFemale;

% Zoom
handles.zoomNr = 0;
handles.radius = radius;
handles.xrange = 1:xmax;
handles.yrange = 1:ymax;
handles.zrange = 1:zmax;

handles.viewNr = 1; % 1: Coronal, 2: Sagittal, 3: Axial
handles.allowDPCB = 1; % allow DragPoint CallBacks?

% Custom colormap for Jacobian
% blue, white, red
handles.BWRmap = BWRmap;


% Start with Coronal on all axes

% axes 1 : undeformed
handles.sparefCor = imref2d([ymax xmax]); % spatial reference
handles.slice1 = squeeze( handles.undefvol(:,:,handles.zval) );
handles.hImage1 = imshow(handles.slice1', handles.sparefCor, 'Parent', handles.axes1 , 'DisplayRange', [0 1]);
set(handles.axes1, 'DataAspectRatio', corscale);
set(handles.hImage1, 'ButtonDownFcn', @axes1_ButtonDownFcn); % set Button Down
handles.colbar1 = colorbar(handles.axes1, 'westOutside');
handles.axes1.XTick = [];
handles.axes1.YTick = [];
handles.dragpoint1 = new_dragpointCor(handles, handles.axes1); % dragpoint

% axes 2 : deformed
handles.slice2 = squeeze( handles.deformvol(:,:,handles.zval) );
handles.hImage2 = imshow(handles.slice2', handles.sparefCor, 'Parent', handles.axes2  , 'DisplayRange' , [0 1]);
set(handles.axes2, 'DataAspectRatio', corscale);
set(handles.hImage2, 'ButtonDownFcn', @axes2_ButtonDownFcn); % set Button Down
handles.colbar2 = colorbar(handles.axes2, 'westOutside');
handles.axes2.XTick = [];
handles.axes2.YTick = [];
handles.dragpoint2 = new_dragpointCor(handles, handles.axes2); % dragpoint

% axes 3 : Jacobian Determinant
handles.slice3 = squeeze( handles.jacdetvol(:,:,handles.zval) );
handles.hImage3 = imshow(handles.slice3', handles.sparefCor, 'Parent', handles.axes3 , 'DisplayRange' , [0 1]);
set(handles.axes3, 'DataAspectRatio', corscale);
set(handles.hImage3, 'ButtonDownFcn', @axes3_ButtonDownFcn); % set Button Down
handles.colbar3 = colorbar(handles.axes3, 'eastOutside');
handles.axes3.XTick = [];
handles.axes3.YTick = [];
handles.dragpoint3 = new_dragpointCor(handles, handles.axes3); % dragpoint
colormap(handles.axes3, handles.BWRmap );
set(handles.axes3, 'CLim', [0 2]);

% axes 4 : Reference volume
handles.slice4 = squeeze( handles.refvolFemale(:,:,handles.zval) );
handles.hImage4 = imshow(handles.slice4', handles.sparefCor, 'Parent', handles.axes4 , 'DisplayRange' , [0 1]);
set(handles.axes4, 'DataAspectRatio', corscale);
set(handles.hImage4, 'ButtonDownFcn', @axes4_ButtonDownFcn); % set Button Down
handles.colbar4 = colorbar(handles.axes4, 'westOutside');
handles.axes4.XTick = [];
handles.axes4.YTick = [];
handles.dragpoint4 = new_dragpointCor(handles, handles.axes4); % dragpoint

% axes 5 : Difference volume
handles.slice5 = squeeze( handles.diffvol(:,:,handles.zval) );
handles.hImage5 = imshow(handles.slice5', handles.sparefCor, 'Parent', handles.axes5 , 'DisplayRange', [0 1]);
set(handles.axes5, 'DataAspectRatio', corscale);
set(handles.hImage5, 'ButtonDownFcn', @axes5_ButtonDownFcn); % set Button Down
handles.colbar5 = colorbar(handles.axes5, 'westOutside');
handles.axes5.XTick = [];
handles.axes5.YTick = [];
handles.dragpoint5 = new_dragpointCor(handles, handles.axes5); % dragpoint
colormap(handles.axes5, handles.BWRmap );
set(handles.axes5, 'CLim', [-2 2]);

% Choose default command line output for del2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes del2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = del2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function editID_Callback(hObject, eventdata, handles)
% hObject    handle to editID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
new_value = round(str2double(get(hObject, 'String')));
handles.poemID = new_value;
guidata(hObject, handles);
% update axes
load_volumes(handles);
update_axes();
% Hints: get(hObject,'String') returns contents of editID as text
%        str2double(get(hObject,'String')) returns contents of editID as a double


% --- Executes during object creation, after setting all properties.
function editID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonLoad.
function pushbuttonLoad_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.test
    % no change
else
load('IDselected.mat');
handles.poemID = IDselected;
set(handles.editID, 'String', int2str(handles.poemID) );
guidata(hObject, handles);
% update axes
load_volumes(handles);
update_axes();
end


function load_volumes(handles)
IDstring = int2str(handles.poemID);

D = dir(handles.pathDeformed);
U = dir(handles.pathUndeformed);
sizeD = length(D);
sizeU = length(U);
deformNotFound = 1;
jacdetNotFound = 1;
undefNotFound = 1;

for i = 1 : sizeD
    strr = D(i).name;
    if contains(strr, 'fat_') && contains(strr, strcat(IDstring, '.vtk'))
        info = vtk_read_header(strr);
        deformvol = vtk_read_volume(info);
        deformNotFound = 0;
    end
    if contains(strr, 'JacDet_') && contains(strr, strcat(IDstring, '.vtk'))
        info = vtk_read_header(strr);
        jacdetvol = vtk_read_volume(info);
        jacdetNotFound = 0;
    end
end

for i = 1 : sizeU
    strr = U(i).name;
    if contains(strr, strcat(IDstring, '_fat_'))
        % read volume
        info = vtk_read_header(strr);
        undefvol = vtk_read_volume(info);
        undefNotFound = 0;
    end
end

if deformNotFound
    sprintf(strcat("No deformed volume found for ID ", IDstring))
    deformvol = zeros(handles.xmax, handles.ymax, handles.zmax);
end
    
if jacdetNotFound
    sprintf(strcat("No Jacobian Determinant volume found for ID ", IDstring))
    jacdetvol = zeros(handles.xmax, handles.ymax, handles.zmax);
end

if undefNotFound
    sprintf(strcat("No undeformed volume found for ID ", IDstring))
    undefvol = zeros(handles.xmax, handles.ymax, handles.zmax);
end

handles.deformvol = deformvol;
handles.jacdetvol = jacdetvol;
handles.undefvol = undefvol;

i = 1;
while i <= handles.poemlength
    if handles.poemID == handles.poemdxa{i, 1}
        handles.female = handles.poemdxa{i, 2};
        i = handles.poemlength + 1;
    else
        i = i + 1;
    end
end

if handles.female
    handles.diffvol = deformvol - handles.refvolFemale;
else
    handles.diffvol = deformvol - handles.refvolMale;
end

%sprintf("load_volumes done")

guidata(gcf, handles);


function update_axes()

handles = guidata(gcf);
% Needed if zoomed
xrange = handles.xrange;
yrange = handles.yrange;
zrange = handles.zrange;

if handles.female
    currentrefvol = handles.refvolFemale;
    set(handles.textRef, 'String', 'Reference Subject: Female');
    set(handles.textRef_2, 'String', 'Reference Subject: Female');
else
    currentrefvol = handles.refvolMale;
    set(handles.textRef, 'String', 'Reference Subject: Male');
    set(handles.textRef_2, 'String', 'Reference Subject: Male');
end


switch handles.viewNr
    case 1
        % Coronal
        handles.slice1 = squeeze( handles.undefvol(xrange,yrange,handles.zval) );
        handles.slice2 = squeeze( handles.deformvol(xrange,yrange,handles.zval) );
        handles.slice3 = squeeze( handles.jacdetvol(xrange,yrange,handles.zval) );
        handles.slice4 = squeeze( currentrefvol(xrange,yrange,handles.zval) );
        handles.slice5 = squeeze( handles.diffvol(xrange,yrange,handles.zval) );
        set(handles.hImage1, 'CData', handles.slice1');
        set(handles.hImage2, 'CData', handles.slice2');
        set(handles.hImage3, 'CData', handles.slice3');
        set(handles.hImage4, 'CData', handles.slice4');
        set(handles.hImage5, 'CData', handles.slice5');
        set(handles.axes1, 'DataAspectRatio', handles.corscale);
        set(handles.axes2, 'DataAspectRatio', handles.corscale);
        set(handles.axes3, 'DataAspectRatio', handles.corscale);
        set(handles.axes4, 'DataAspectRatio', handles.corscale);
        set(handles.axes5, 'DataAspectRatio', handles.corscale);
        guidata(gcf, handles);
        update_dragpointsCor(handles);
    case 2
        % Sagittal
        handles.slice1 = squeeze( handles.undefvol(handles.xval,yrange,zrange) );
        handles.slice2 = squeeze( handles.deformvol(handles.xval,yrange,zrange) );
        handles.slice3 = squeeze( handles.jacdetvol(handles.xval,yrange,zrange) );
        handles.slice4 = squeeze( currentrefvol(handles.xval,yrange,zrange) );
        handles.slice5 = squeeze( handles.diffvol(handles.xval,yrange,zrange) );
        set(handles.hImage1, 'CData', handles.slice1);
        set(handles.hImage2, 'CData', handles.slice2);
        set(handles.hImage3, 'CData', handles.slice3);
        set(handles.hImage4, 'CData', handles.slice4);
        set(handles.hImage5, 'CData', handles.slice5);
        set(handles.axes1, 'DataAspectRatio', handles.sagscale);
        set(handles.axes2, 'DataAspectRatio', handles.sagscale);
        set(handles.axes3, 'DataAspectRatio', handles.sagscale);
        set(handles.axes4, 'DataAspectRatio', handles.sagscale);
        set(handles.axes5, 'DataAspectRatio', handles.sagscale);
        guidata(gcf, handles);
        update_dragpointsSag(handles);
    case 3
        % Axial
        handles.slice1 = squeeze( handles.undefvol(xrange,handles.yval,zrange) );
        handles.slice2 = squeeze( handles.deformvol(xrange,handles.yval,zrange) );
        handles.slice3 = squeeze( handles.jacdetvol(xrange,handles.yval,zrange) );
        handles.slice4 = squeeze( currentrefvol(xrange,handles.yval,zrange) );
        handles.slice5 = squeeze( handles.diffvol(xrange,handles.yval,zrange) );
        set(handles.hImage1, 'CData', handles.slice1');
        set(handles.hImage2, 'CData', handles.slice2');
        set(handles.hImage3, 'CData', handles.slice3');
        set(handles.hImage4, 'CData', handles.slice4');
        set(handles.hImage5, 'CData', handles.slice5');
        set(handles.axes1, 'DataAspectRatio', handles.axiscale);
        set(handles.axes2, 'DataAspectRatio', handles.axiscale);
        set(handles.axes3, 'DataAspectRatio', handles.axiscale);
        set(handles.axes4, 'DataAspectRatio', handles.axiscale);
        set(handles.axes5, 'DataAspectRatio', handles.axiscale);
        guidata(gcf, handles);
        update_dragpointsAxi(handles);
end
%sprintf("update_axes done")

update_values();
update_xyz();
update_slider();



% --- Executes on button press in radiobuttonCor.
function radiobuttonCor_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonCor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.viewNr = 1;
guidata(gcf, handles);
update_axes();
% Hint: get(hObject,'Value') returns toggle state of radiobuttonCor


% --- Executes on button press in radiobuttonSag.
function radiobuttonSag_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonSag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.viewNr = 2;
guidata(gcf, handles);
update_axes();
% Hint: get(hObject,'Value') returns toggle state of radiobuttonSag


% --- Executes on button press in radiobuttonAxi.
function radiobuttonAxi_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonAxi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.viewNr = 3;
guidata(gcf, handles);
update_axes();
% Hint: get(hObject,'Value') returns toggle state of radiobuttonAxi


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
            z = handles.zval + 1;
            handles.zval = coordlimit(z, handles.zmax);
            guidata(hObject, handles);
            update_axes();
        case 2
            x = handles.xval + 1;
            handles.xval = coordlimit(x, handles.xmax);
            guidata(hObject, handles);
            update_axes();
        case 3
            y = handles.yval + 1;
            handles.yval = coordlimit(y, handles.ymax);
            guidata(hObject, handles);
            update_axes();
    end
else
    switch handles.viewNr
        case 1
            z = handles.zval - 1;
            handles.zval = coordlimit(z, handles.zmax);
            guidata(hObject, handles);
            update_axes();
        case 2
            x = handles.xval - 1;
            handles.xval = coordlimit(x, handles.xmax); 
            guidata(hObject, handles);
            update_axes();
        case 3
            y = handles.yval - 1;
            handles.yval = coordlimit(y, handles.ymax);
            guidata(hObject, handles);
            update_axes();
    end
end


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
    case 'c'
        radiobuttonCor_Callback(hObject, eventdata, handles);
        set(handles.radiobuttonCor, 'Value', 1);
    case 's'
        radiobuttonSag_Callback(hObject, eventdata, handles);
        set(handles.radiobuttonSag, 'Value', 1);
    case 'a'
        radiobuttonAxi_Callback(hObject, eventdata, handles);
        set(handles.radiobuttonAxi, 'Value', 1);
    case 'l'
        uicontrol(handles.pushbuttonLoad);
        pushbuttonLoad_Callback(hObject, eventdata, handles);
    case 'i'
        uicontrol(handles.pushbuttonZoomIn);
        pushbuttonZoomIn_Callback(hObject, eventdata, handles);
    case 'o'
        uicontrol(handles.pushbuttonZoomOut);
        pushbuttonZoomOut_Callback(hObject, eventdata, handles);
    case 'x'
        uicontrol(handles.editX);
    case 'y'
        uicontrol(handles.editY);
    case 'z'
        uicontrol(handles.editZ);   
    case 'f'
        uicontrol(handles.pushbuttonFindLimits);
        pushbuttonFindLimits_Callback(hObject, eventdata, handles);
        
    otherwise
        %sprintf("That key is unbound.")
end


% Make New and Update Dragpoints

function dp = new_dragpointCor(handles, axesObj)
dp = impoint(axesObj, handles.xval, handles.yval);
setColor(dp, 'magenta'); %black
dpCor = @dragpointCallback;
addNewPositionCallback(dp, dpCor);

function update_dragpointsCor(handles)
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
setPosition(handles.dragpoint1, xval, yval);
setPosition(handles.dragpoint2, xval, yval);
setPosition(handles.dragpoint3, xval, yval);
setPosition(handles.dragpoint4, xval, yval);
setPosition(handles.dragpoint5, xval, yval);
handles.allowDPCB = 1;
guidata(gcf, handles);

function update_dragpointsSag(handles)
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
setPosition(handles.dragpoint1, zval, yval);
setPosition(handles.dragpoint2, zval, yval);
setPosition(handles.dragpoint3, zval, yval);
setPosition(handles.dragpoint4, zval, yval);
setPosition(handles.dragpoint5, zval, yval);
handles.allowDPCB = 1;
guidata(gcf, handles);

function update_dragpointsAxi(handles)
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
setPosition(handles.dragpoint1, xval, zval);
setPosition(handles.dragpoint2, xval, zval);
setPosition(handles.dragpoint3, xval, zval);
setPosition(handles.dragpoint4, xval, zval);
setPosition(handles.dragpoint5, xval, zval);
handles.allowDPCB = 1;
guidata(gcf, handles);


% Dragpoint Callbacks

% --- Exectes when dragpoint is moved (in ANY axis 1-4)
% pos = [x y]
function dragpointCallback(pos)
handles = guidata(gcf);
znr = handles.zoomNr;
xrange = handles.xrange;
yrange = handles.yrange;
zrange = handles.zrange;
if handles.allowDPCB

    switch handles.viewNr
        case 1 % Cor
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
            update_dragpointsCor(handles);    
        case 2 % Sag
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
            update_dragpointsSag(handles);
        case 3 % Axi
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
            update_dragpointsAxi(handles);
    end
    guidata(gcf, handles);
    
    update_xyz();
    update_values();
    
else
    % no action
    %sprintf("Callback Cor disallowed")
end


function update_values()
handles = guidata(gcf);
x = handles.xval;
y = handles.yval;
z = handles.zval;

string1 = num2str(handles.undefvol(x,y,z));
string2 = num2str(handles.deformvol(x,y,z));
string3 = num2str(handles.jacdetvol(x,y,z));

set(handles.textVal1_2, 'String', string1);
set(handles.textVal2_2, 'String', string2);
set(handles.textVal3_2, 'String', string3);

if handles.female
    currentrefvol = handles.refvolFemale;
else
    currentrefvol = handles.refvolMale;
end
string4 = num2str(currentrefvol(x,y,z));
set(handles.textVal4_2, 'String', string4);

guidata(gcf, handles);




function editX_Callback(hObject, eventdata, handles)
% hObject    handle to editX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
new_xval = round(str2double(get(hObject, 'String')));

if new_xval > 0 && new_xval <= handles.xmax
    handles.xval = new_xval;
    
    guidata(hObject, handles);
    update_axes();
    
    switch handles.viewNr
        case 1 % Cor
            update_dragpointsCor(handles);    
        case 2 % Sag
            update_dragpointsSag(handles);
        case 3 % Axi
            update_dragpointsAxi(handles);
    end

else
    % no action
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



function editY_Callback(hObject, eventdata, handles)
% hObject    handle to editY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
new_yval = round(str2double(get(hObject, 'String')));

if new_yval > 0 && new_yval <= handles.ymax
    handles.yval = new_yval;
    
    guidata(hObject, handles);
    update_axes();
    
    switch handles.viewNr
        case 1 % Cor
            update_dragpointsCor(handles);    
        case 2 % Sag
            update_dragpointsSag(handles);
        case 3 % Axi
            update_dragpointsAxi(handles);
    end
    guidata(hObject, handles);
else
    % no action
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



function editZ_Callback(hObject, eventdata, handles)
% hObject    handle to editZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
new_zval = round(str2double(get(hObject, 'String')));

if new_zval > 0 && new_zval <= handles.zmax
    handles.zval = new_zval;

    guidata(hObject, handles);
    update_axes();
    
    switch handles.viewNr
        case 1 % Cor
            update_dragpointsCor(handles);    
        case 2 % Sag
            update_dragpointsSag(handles);
        case 3 % Axi
            update_dragpointsAxi(handles);
    end
    guidata(hObject, handles);
else
    % no action
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


function update_xyz()
handles = guidata(gcf);
x = handles.xval;
y = handles.yval;
z = handles.zval;

set(handles.editX, 'String', int2str(x));
set(handles.editY, 'String', int2str(y));
set(handles.editZ, 'String', int2str(z));

guidata(gcf, handles);




% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%min get(hObject, 'Min') = 0
%max get(hObject, 'Max') = 1
value = get(hObject, 'Value');

%handles.viewNr = 1: Coronal(z), 2: Sagittal(x), 3: Axial(y)
switch handles.viewNr
    case 1
        new_z = max(round(handles.zmax * value), 1); % z min = 1
        handles.zval = new_z;
    case 2
        new_x = max(round(handles.xmax * value), 1); % x min = 1
        handles.xval = new_x;
    case 3
        new_y = max(round(handles.ymax * value), 1); % y min = 1
        handles.yval = new_y;
end
guidata(hObject, handles);
update_axes();
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function update_slider()
handles = guidata(gcf);

switch handles.viewNr
    case 1
        new_z = handles.zval; 
        value = (new_z - 1) / handles.zmax; % z min = 1 
        set(handles.slider1, 'Value', value);
    case 2
        new_x = handles.xval; 
        value = (new_x - 1) / handles.xmax; % x min = 1 
        set(handles.slider1, 'Value', value);
    case 3
        new_y = handles.yval;
        value = (new_y - 1) / handles.ymax; % y min = 1 
        set(handles.slider1, 'Value', value);
end
guidata(gcf, handles);


% --- Executes on button press in pushbuttonReset.
function pushbuttonReset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Reset all coordinates to middle
handles.xval = handles.xinit;
handles.yval = handles.yinit;
handles.zval = handles.zinit;

% Reset Zoom
handles.xrange = 1:handles.xmax;
handles.yrange = 1:handles.ymax;
handles.zrange = 1:handles.zmax;
handles.zoomNr = 0;

guidata(hObject, handles);
update_axes();
switch handles.viewNr
    case 1
        update_dragpointsCor(handles);
    case 2
        update_dragpointsSag(handles);
    case 3
        update_dragpointsAxi(handles);
end


% JD limits

function editJDupper_Callback(hObject, eventdata, handles)
% hObject    handle to editJDupper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = str2double(get(hObject, 'String'));
lowerupper = get(handles.axes3, 'CLim');

if value > 5 % 10
    value = 5;
end
if value < 1
    value = 1;
end
set(handles.axes3, 'CLim', [lowerupper(1) value]);
set(handles.sliderJDupper, 'Value', (value - 1)/4); % /10

guidata(hObject, handles);
update_axes();
% Hints: get(hObject,'String') returns contents of editJDupper as text
%        str2double(get(hObject,'String')) returns contents of editJDupper as a double


% --- Executes during object creation, after setting all properties.
function editJDupper_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editJDupper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderJDupper_Callback(hObject, eventdata, handles)
% hObject    handle to sliderJDupper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lowerupper = get(handles.axes3, 'CLim');
% value [0 1]
value = get(hObject, 'Value');
new_upper = value * 4 + 1; % * 10
set(handles.axes3, 'CLim', [lowerupper(1) new_upper]);
set(handles.editJDupper, 'String', num2str(new_upper));
guidata(hObject, handles);
update_axes();
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sliderJDupper_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderJDupper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function editJDlower_Callback(hObject, eventdata, handles)
% hObject    handle to editJDlower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = str2double(get(hObject, 'String'));
lowerupper = get(handles.axes3, 'CLim');

if value < -3 % 10
    value = -3;
end
if value >= 1
    value = 0.99;
end
set(handles.axes3, 'CLim', [value lowerupper(2)]);
set(handles.sliderJDlower, 'Value', (value + 3)/4); % + 10/10
guidata(hObject, handles);
update_axes();
% Hints: get(hObject,'String') returns contents of editJDlower as text
%        str2double(get(hObject,'String')) returns contents of editJDlower as a double


% --- Executes during object creation, after setting all properties.
function editJDlower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editJDlower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderJDlower_Callback(hObject, eventdata, handles)
% hObject    handle to sliderJDlower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lowerupper = get(handles.axes3, 'CLim');
% value [0 1]
value = get(hObject, 'Value');
new_lower = (value - 1) * 4 + 1;

if new_lower > 0.99
    new_lower = 0.99;
end

set(handles.axes3, 'CLim', [new_lower lowerupper(2)]);
set(handles.editJDlower, 'String', num2str(new_lower));
guidata(hObject, handles);
update_axes();
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sliderJDlower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderJDlower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% Limits x,y,z coords from being < 1 or > max 
function aout = coordlimit(a, alim)
aout = a;
if a < 1
    aout = 1;
end
if a > alim
    aout = alim;
end


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
    update_axes();
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
    update_axes();
end


% Popup Menus

% --- Executes on selection change in popupmenuColormap.
function popupmenuColormap_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuColormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(hObject,'Value');
switch value
    case 1
        colormap(handles.axes1, jet );
        colormap(handles.axes2, jet );
        colormap(handles.axes4, jet );
    case 2
        colormap(handles.axes1, gray );
        colormap(handles.axes2, gray );
        colormap(handles.axes4, gray );
    case 3
        colormap(handles.axes1, parula );
        colormap(handles.axes2, parula );
        colormap(handles.axes4, parula );
    case 4
        colormap(handles.axes1, hot );
        colormap(handles.axes2, hot );
        colormap(handles.axes4, hot );
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


% --- Executes on selection change in popupmenuMarkerColor.
function popupmenuMarkerColor_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuMarkerColor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(hObject, 'Value'); 
switch value
    case 1
    setColor(handles.dragpoint1, 'black');
    setColor(handles.dragpoint2, 'black');
    setColor(handles.dragpoint3, 'black');
    setColor(handles.dragpoint4, 'black');
    setColor(handles.dragpoint5, 'black');
    case 2
    setColor(handles.dragpoint1, 'white');
    setColor(handles.dragpoint2, 'white');
    setColor(handles.dragpoint3, 'white');
    setColor(handles.dragpoint4, 'white');
    setColor(handles.dragpoint5, 'white');
    case 3
    setColor(handles.dragpoint1, 'magenta');
    setColor(handles.dragpoint2, 'magenta');
    setColor(handles.dragpoint3, 'magenta');
    setColor(handles.dragpoint4, 'magenta');
    setColor(handles.dragpoint5, 'magenta');
end
guidata(gcf, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuMarkerColor contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuMarkerColor


% --- Executes during object creation, after setting all properties.
function popupmenuMarkerColor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuMarkerColor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenuColormapDiff.
function popupmenuColormapDiff_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuColormapDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(hObject,'Value');
switch value
    case 1 % BWR
        colormap(handles.axes5, handles.BWRmap );
    case 2 % Jet
        colormap(handles.axes5, jet );
    case 3 % Gray
        colormap(handles.axes5, gray );
end
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuColormapDiff contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuColormapDiff


% --- Executes during object creation, after setting all properties.
function popupmenuColormapDiff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuColormapDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuColormapJD.
function popupmenuColormapJD_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuColormapJD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = get(hObject,'Value');
switch value
    case 1 % BWR
        colormap(handles.axes3, handles.BWRmap );
    case 2 % Jet
        colormap(handles.axes3, jet );
    case 3 % Gray
        colormap(handles.axes3, gray );
end
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuColormapJD contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuColormapJD


% --- Executes during object creation, after setting all properties.
function popupmenuColormapJD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuColormapJD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% Button Down Functions : Handle clicking on images

% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

point = get(handles.axes1, 'CurrentPoint');
imX = point(1, 1);
imY = point(1, 2);

% The unsupressed callback function takes care of updating everything
setPosition(handles.dragpoint1, imX, imY);
setPosition(handles.dragpoint2, imX, imY);
setPosition(handles.dragpoint3, imX, imY);
setPosition(handles.dragpoint4, imX, imY);
setPosition(handles.dragpoint5, imX, imY);


% --- Executes on mouse press over axes background.
function axes2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

point = get(handles.axes2, 'CurrentPoint');
imX = point(1, 1);
imY = point(1, 2);

setPosition(handles.dragpoint1, imX, imY);
setPosition(handles.dragpoint2, imX, imY);
setPosition(handles.dragpoint3, imX, imY);
setPosition(handles.dragpoint4, imX, imY);
setPosition(handles.dragpoint5, imX, imY);


% --- Executes on mouse press over axes background.
function axes3_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

point = get(handles.axes3, 'CurrentPoint');
imX = point(1, 1);
imY = point(1, 2);

setPosition(handles.dragpoint1, imX, imY);
setPosition(handles.dragpoint2, imX, imY);
setPosition(handles.dragpoint3, imX, imY);
setPosition(handles.dragpoint4, imX, imY);
setPosition(handles.dragpoint5, imX, imY);


% --- Executes on mouse press over axes background.
function axes4_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

point = get(handles.axes4, 'CurrentPoint');
imX = point(1, 1);
imY = point(1, 2);

setPosition(handles.dragpoint1, imX, imY);
setPosition(handles.dragpoint2, imX, imY);
setPosition(handles.dragpoint3, imX, imY);
setPosition(handles.dragpoint4, imX, imY);
setPosition(handles.dragpoint5, imX, imY);


% --- Executes on mouse press over axes background.
function axes5_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

point = get(handles.axes5, 'CurrentPoint');
imX = point(1, 1);
imY = point(1, 2);

setPosition(handles.dragpoint1, imX, imY);
setPosition(handles.dragpoint2, imX, imY);
setPosition(handles.dragpoint3, imX, imY);
setPosition(handles.dragpoint4, imX, imY);
setPosition(handles.dragpoint5, imX, imY);



% --- Executes on button press in pushbuttonFindLimits.
function pushbuttonFindLimits_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonFindLimits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
maxval = max(max(handles.slice3));
minval = min(min(handles.slice3));

lower = minval;
upper = maxval;

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

set(handles.editJDlower, 'String', num2str(lower));
set(handles.editJDupper, 'String', num2str(upper));
set(handles.sliderJDlower, 'Value', (lower + 3) / 4);
set(handles.sliderJDupper, 'Value', (upper - 1) / 4);

set(handles.axes3, 'CLim', [lower upper]);


% --- Executes on button press in pushbuttonResetLims.
function pushbuttonResetLims_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonResetLims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lower = 0;
upper = 2;

set(handles.editJDlower, 'String', num2str(lower));
set(handles.editJDupper, 'String', num2str(upper));
set(handles.sliderJDlower, 'Value', (lower + 3) / 4);
set(handles.sliderJDupper, 'Value', (upper - 1) / 4);

set(handles.axes3, 'CLim', [lower upper]);

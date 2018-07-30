function varargout = Simulator7vol2(varargin)

% Last Modified by GUIDE v2.5 18-May-2018 18:42:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Simulator7vol2_OpeningFcn, ...
                   'gui_OutputFcn',  @Simulator7vol2_OutputFcn, ...
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



% --- Executes just before Simulator7vol2 is made visible.
function Simulator7vol2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)


axes(handles.axes4);
I = imread('AIT_full.png');
imshow(I);

axes(handles.axes5);
J = imread('B_WiSE_logo_1full_white.png');
imshow(J);

axes(handles.axes6);
N = imread('TCD_LOGO.png');
imshow(N);

handles.output = hObject;   % Choose default command line output for Simulator6
guidata(hObject, handles);  % Update handles structure


% --- Outputs from this function are returned to the command line.
function varargout = Simulator7vol2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function simmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to simmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function misprob2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to misprob2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles,data)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
%penaltyThreshold = getappdata(0,'penaltyThreshold')  %get data for every variable from GUI
n_sim = getappdata(0,'n_sim');

Mno11 = getappdata(0,'Mno11');
Mno21 = getappdata(0,'Mno21');

threspen21 = getappdata(0,'threspen21');
threspen22 = getappdata(0,'threspen22');
threspen23 = getappdata(0,'threspen23');


simulatedTime=getappdata(0,'simulatedTime');   %different file than the previous
getInputParametersPenThresh();
fun_runMultipleSimulationsPenThresh(n_sim,Mno11,Mno21,threspen21,threspen22,threspen23);%,threspen24,threspen25);


function numsim2_Callback(hObject, eventdata, handles)
% hObject    handle to numsim2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numsim2 as text
%        str2double(get(hObject,'String')) returns contents of numsim2 as a double

n_sim = str2num(get(handles.numsim2, 'String'));
setappdata(0,'n_sim',n_sim);


% --- Executes during object creation, after setting all properties.
function numsim2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numsim2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in simtime2.
function simtime2_Callback(hObject, eventdata, handles)
% hObject    handle to simtime2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns simtime2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from simtime2

simulatedTime = str2num(get(handles.simtime2, 'String'));   %take simulatedTime from GUI
setappdata(0,'simulatedTime',simulatedTime);


% --- Executes during object creation, after setting all properties.
function simtime2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to simtime2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Mno11_Callback(hObject, eventdata, handles)
% hObject    handle to Mno11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Mno11 as text
%        str2double(get(hObject,'String')) returns contents of Mno11 as a double

Mno11 = str2num(get(handles.Mno11, 'String'));   %take first misbehavior probability for MNO from GUI
setappdata(0,'Mno11',Mno11);


% --- Executes during object creation, after setting all properties.
function Mno11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mno11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Mno21_Callback(hObject, eventdata, handles)
% hObject    handle to Mno21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Mno21 as text
%        str2double(get(hObject,'String')) returns contents of Mno21 as a double

Mno21 = str2num(get(handles.Mno21, 'String'));   %take first misbehavior probability for MNO from GUI
setappdata(0,'Mno21',Mno21);


% --- Executes during object creation, after setting all properties.
function Mno21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mno21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function threspen21_Callback(hObject, eventdata, handles)
% hObject    handle to threspen21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threspen21 as text
%        str2double(get(hObject,'String')) returns contents of threspen21 as a double

threspen21 = str2num(get(handles.threspen21, 'String'));   %take PenaltyThreshold from GUI
setappdata(0,'threspen21',threspen21);


% --- Executes during object creation, after setting all properties.
function threspen21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threspen21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function threspen22_Callback(hObject, eventdata, handles)
% hObject    handle to threspen22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threspen22 as text
%        str2double(get(hObject,'String')) returns contents of threspen22 as a double

threspen22 = str2num(get(handles.threspen22, 'String'));   %take PenaltyThreshold from GUI
setappdata(0,'threspen22',threspen22);


% --- Executes during object creation, after setting all properties.
function threspen22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threspen22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function threspen23_Callback(hObject, eventdata, handles)
% hObject    handle to threspen23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threspen23 as text
%        str2double(get(hObject,'String')) returns contents of threspen23 as a double

threspen23 = str2num(get(handles.threspen23, 'String'));   %take PenaltyThreshold from GUI
setappdata(0,'threspen23',threspen23);

% --- Executes during object creation, after setting all properties.
function threspen23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threspen23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function threspen24_Callback(hObject, eventdata, handles)
% hObject    handle to threspen24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threspen24 as text
%        str2double(get(hObject,'String')) returns contents of threspen24 as a double

threspen24 = str2num(get(handles.threspen24, 'String'));   %take PenaltyThreshold from GUI
setappdata(0,'threspen24',threspen24);

% --- Executes during object creation, after setting all properties.
function threspen24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threspen24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function threspen25_Callback(hObject, eventdata, handles)
% hObject    handle to threspen25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threspen25 as text
%        str2double(get(hObject,'String')) returns contents of threspen25 as a double


threspen25 = str2num(get(handles.threspen25, 'String'));   %take PenaltyThreshold from GUI
setappdata(0,'threspen25',threspen25);

% --- Executes during object creation, after setting all properties.
function threspen25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threspen25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

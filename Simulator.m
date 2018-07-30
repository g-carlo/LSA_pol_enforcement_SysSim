function varargout = Simulator(varargin)
%SIMULATOR MATLAB code file for Simulator.fig
%      SIMULATOR, by itself, creates a new SIMULATOR or raises the existing
%      singleton*.
%
%      H = SIMULATOR returns the handle to a new SIMULATOR or the handle to
%      the existing singleton*.
%
%      SIMULATOR('Property','Value',...) creates a new SIMULATOR using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to Simulator_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      SIMULATOR('CALLBACK') and SIMULATOR('CALLBACK',hObject,...) call the
%      local function named CALLBACK in SIMULATOR.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Simulator

% Last Modified by GUIDE v2.5 20-Apr-2018 09:36:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Simulator_OpeningFcn, ...
                   'gui_OutputFcn',  @Simulator_OutputFcn, ...
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


% --- Executes just before Simulator is made visible.
function Simulator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)
I = imread('C:\Users\bwise-Alpha\Desktop\EFI\1. AIT\LOGOS\logo.png');
imshow(I);
% Choose default command line output for Simulator
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Simulator wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Simulator_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in simmenu2.
function simmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to simmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
g
% Hints: contents = cellstr(get(hObject,'String')) returns simmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from simmenu2


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


% --- Executes on selection change in incumbmenu3.
%function incumbmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to incumbmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns incumbmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from incumbmenu3


% --- Executes during object creation, after setting all properties.
%function incumbmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to incumbmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in misprob2.
function misprob2_Callback(hObject, eventdata, handles)
% hObject    handle to misprob2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns misprob2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from misprob2


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
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% pare tis parametrous apo to gui

numsimulation = str2num(get(handles.numsim2, 'string'));
threspenalty = str2num(get(handles.threspen2, 'string'));
simtime = str2num(get(handles.simtime2, 'string'));



sum = numsimulation + threspenalty + simtime;
s = num2str(sum)

% % if sinthiki gia na doume an tha paroume figure6 h figure7
% 	if
%     runMultipleSimulationsMisbProb.m
%     else
%     runMultipleSimulationsPenThresh.m
%     end
    
% kane tous ipologismous
%sum = n1+n2;

% kane plot ta apotelesmata

%s = num2str(sum);
%set(handles.result,'string',s);


% --- Executes on button press in realtime.
function realtime_Callback(hObject, eventdata, handles)
% hObject    handle to realtime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of realtime


% --- Executes on button press in offline.
function offline_Callback(hObject, eventdata, handles)
% hObject    handle to offline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of offline



function numsim2_Callback(hObject, eventdata, handles)
% hObject    handle to numsim2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numsim2 as text
%        str2double(get(hObject,'String')) returns contents of numsim2 as a double


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



function threspen2_Callback(hObject, eventdata, handles)
% hObject    handle to threspen2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threspen2 as text
%        str2double(get(hObject,'String')) returns contents of threspen2 as a double


% --- Executes during object creation, after setting all properties.
function threspen2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threspen2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

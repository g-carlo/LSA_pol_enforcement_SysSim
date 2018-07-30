function varargout = Simulator6(varargin)

% Last Modified by GUIDE v2.5 09-May-2018 18:57:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Simulator6_OpeningFcn, ...
                   'gui_OutputFcn',  @Simulator6_OutputFcn, ...
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



% --- Executes just before Simulator6 is made visible.
function Simulator6_OpeningFcn(hObject, eventdata, handles, varargin)
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
function varargout = Simulator6_OutputFcn(hObject, eventdata, handles)
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
penaltyGiven = getappdata(0,'penaltyGiven');  %get data for every variable from GUI
n_sim = getappdata(0,'n_sim');

Mno11 = getappdata(0,'Mno11');
Mno12 = getappdata(0,'Mno12');
Mno13 = getappdata(0,'Mno13');
Mno14 = getappdata(0,'Mno14');

Mno21 = getappdata(0,'Mno21');
Mno22 = getappdata(0,'Mno22');
Mno23 = getappdata(0,'Mno23');
Mno24 = getappdata(0,'Mno24');

%simulatedTime_gui=getappdata(0,'simulatedTime');   %different file than the previous

getInputParametersMisProb();
fun_runMultipleSimulationsMisProb(penaltyGiven,n_sim,Mno11,Mno12,Mno13,Mno14,Mno21,Mno22,Mno23,Mno24);

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


% % --- Executes on selection change in simtime2.
% function simtime2_Callback(hObject, eventdata, handles)
% % hObject    handle to simtime2 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: contents = cellstr(get(hObject,'String')) returns simtime2 contents as cell array
% %        contents{get(hObject,'Value')} returns selected item from simtime2
% 
% simulatedTime = str2num(get(handles.simtime2, 'String'));   %take simulatedTime from GUI
% setappdata(0,'simulatedTime',simulatedTime);
% 
% 
% % --- Executes during object creation, after setting all properties.
% function simtime2_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to simtime2 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: popupmenu controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end


function threspen2_Callback(hObject, eventdata, handles)
% hObject    handle to threspen2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threspen2 as text
%        str2double(get(hObject,'String')) returns contents of threspen2 as a double

penaltyGiven = str2num(get(handles.threspen2, 'String'));   %take penaltyGiven from GUI
setappdata(0,'penaltyGiven',penaltyGiven);


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



function Mno12_Callback(hObject, eventdata, handles)
% hObject    handle to Mno12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Mno12 as text
%        str2double(get(hObject,'String')) returns contents of Mno12 as a double

Mno12 = str2num(get(handles.Mno12, 'String'));   %take first misbehavior probability for MNO from GUI
setappdata(0,'Mno12',Mno12);


% --- Executes during object creation, after setting all properties.
function Mno12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mno12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Mno13_Callback(hObject, eventdata, handles)
% hObject    handle to Mno13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Mno13 as text
%        str2double(get(hObject,'String')) returns contents of Mno13 as a double

Mno13 = str2num(get(handles.Mno13, 'String'));   %take first misbehavior probability for MNO from GUI
setappdata(0,'Mno13',Mno13);



% --- Executes during object creation, after setting all properties.
function Mno13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mno13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Mno14_Callback(hObject, eventdata, handles)
% hObject    handle to Mno14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Mno14 as text
%        str2double(get(hObject,'String')) returns contents of Mno14 as a double

Mno14 = str2num(get(handles.Mno14, 'String'));   %take first misbehavior probability for MNO from GUI
setappdata(0,'Mno14',Mno14);



% --- Executes during object creation, after setting all properties.
function Mno14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mno14 (see GCBO)
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



function Mno22_Callback(hObject, eventdata, handles)
% hObject    handle to Mno22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Mno22 as text
%        str2double(get(hObject,'String')) returns contents of Mno22 as a double

Mno22 = str2num(get(handles.Mno22, 'String'));   %take first misbehavior probability for MNO from GUI
setappdata(0,'Mno22',Mno22);


% --- Executes during object creation, after setting all properties.
function Mno22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mno22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Mno23_Callback(hObject, eventdata, handles)
% hObject    handle to Mno23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Mno23 as text
%        str2double(get(hObject,'String')) returns contents of Mno23 as a double

Mno23 = str2num(get(handles.Mno23, 'String'));   %take first misbehavior probability for MNO from GUI
setappdata(0,'Mno23',Mno23);


% --- Executes during object creation, after setting all properties.
function Mno23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mno23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Mno24_Callback(hObject, eventdata, handles)
% hObject    handle to Mno24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Mno24 as text
%        str2double(get(hObject,'String')) returns contents of Mno24 as a double

Mno24 = str2num(get(handles.Mno24, 'String'));   %take first misbehavior probability for MNO from GUI
setappdata(0,'Mno24',Mno24);


% --- Executes during object creation, after setting all properties.
function Mno24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mno24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%get(0, 'ScreenSize');
%set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

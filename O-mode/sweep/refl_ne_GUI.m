                                                                                                                                                                                function varargout = refl_ne_GUI(varargin)
% REFL_NE_GUI M-file for refl_ne_GUI.fig
%      REFL_NE_GUI, by itself, creates a new REFL_NE_GUI or raises the existing
%      singleton*.
%
%      H = REFL_NE_GUI returns the handle to a new REFL_NE_GUI or the handle to
%      the existing singleton*.
%
%      REFL_NE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REFL_NE_GUI.M with the given input arguments.
%
%      REFL_NE_GUI('Property','Value',...) creates a new REFL_NE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before refl_ne_GUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to refl_ne_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help refl_ne_GUI

% Last Modified by GUIDE v2.5 16-Jan-2016 10:13:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @refl_ne_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @refl_ne_GUI_OutputFcn, ...
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


% --- Executes just before refl_ne_GUI is made visible.
function refl_ne_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to refl_ne_GUI (see VARARGIN)

% Choose default command line output for refl_ne_GUI
clc
fid=fopen('refl_init.mat');
if fid==-1
   shot_1='48059';
   IT_1='8000';
   time_1='3';
   data_path_1='W:\2014Exp_ne';
   time_delay_1='20';
   time_trigger_1='-10';
   fs_V_1='1.875';
   fs_IQ_1='60';
   T_period_1='50';
   T_scan_1='40';
   filter_freq_Q_1='[2 8]';
   filter_freq_V_1='[2 10]';
   filter_freq_W_1='[2 12]';
   filter_ICRF_1='25';
   edge_profile_1='1';
   edge_profile_function_1='1';
   judge_cutoff_1='1';
   freq_prob_1='[15:18]';
   rand_number_check_1='1';
   rand_number_1='20';
   judge_Rzero_1='1';
   R_zero_1='2.35';
   save('refl_init.mat','shot_1','IT_1','time_1','data_path_1','time_delay_1', ... 
       'time_trigger_1','fs_V_1','fs_IQ_1','T_period_1','T_scan_1','filter_freq_Q_1','filter_freq_V_1', ...
       'filter_freq_W_1','filter_ICRF_1','edge_profile_1','edge_profile_function_1','judge_cutoff_1',...
   'freq_prob_1','rand_number_check_1','rand_number_1','judge_Rzero_1','R_zero_1');
end
fclose all
load('refl_init.mat');
set(handles.shot,'string',shot_1);
set(handles.IT,'string',IT_1);
set(handles.time,'string',time_1);
set(handles.data_path,'string',data_path_1);
set(handles.time_delay,'string',time_delay_1);
set(handles.time_trigger,'string',time_trigger_1);
set(handles.fs_V,'string',fs_V_1);
set(handles.fs_IQ,'string',fs_IQ_1);
set(handles.T_period,'string',T_period_1);
set(handles.T_scan,'string',T_scan_1);
set(handles.filter_freq_Q,'string',filter_freq_Q_1);
set(handles.filter_freq_V,'string',filter_freq_V_1);
set(handles.filter_freq_W,'string',filter_freq_W_1);
set(handles.filter_ICRF,'string',filter_ICRF_1);
set(handles.edge_profile,'value',str2num(edge_profile_1));
set(handles.edge_profile_function,'value',str2num(edge_profile_function_1));
set(handles.judge_cutoff,'value',str2num(judge_cutoff_1));
set(handles.freq_prob,'string',freq_prob_1);
set(handles.rand_number_check,'value',str2num(rand_number_check_1));
set(handles.rand_number,'string',rand_number_1);
set(handles.judge_Rzero,'value',str2num(judge_Rzero_1));
set(handles.R_zero,'string',R_zero_1);

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes refl_ne_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = refl_ne_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function filter_freq_Q_Callback(hObject, eventdata, handles)
% hObject    handle to filter_freq_Q (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filter_freq_Q as text
%        str2double(get(hObject,'String')) returns contents of filter_freq_Q as a double


% --- Executes during object creation, after setting all properties.
function filter_freq_Q_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filter_freq_Q (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function filter_freq_V_Callback(hObject, eventdata, handles)
% hObject    handle to filter_freq_V (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filter_freq_V as text
%        str2double(get(hObject,'String')) returns contents of filter_freq_V as a double


% --- Executes during object creation, after setting all properties.
function filter_freq_V_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filter_freq_V (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function filter_freq_W_Callback(hObject, eventdata, handles)
% hObject    handle to filter_freq_W (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filter_freq_W as text
%        str2double(get(hObject,'String')) returns contents of filter_freq_W as a double


% --- Executes during object creation, after setting all properties.
function filter_freq_W_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filter_freq_W (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function shot_Callback(hObject, eventdata, handles)
% hObject    handle to shot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of shot as text
%        str2double(get(hObject,'String')) returns contents of shot as a double


% --- Executes during object creation, after setting all properties.
function shot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to shot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function IT_Callback(hObject, eventdata, handles)
% hObject    handle to IT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IT as text
%        str2double(get(hObject,'String')) returns contents of IT as a double


% --- Executes during object creation, after setting all properties.
function IT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function time_Callback(hObject, eventdata, handles)
% hObject    handle to time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time as text
%        str2double(get(hObject,'String')) returns contents of time as a double


% --- Executes during object creation, after setting all properties.
function time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function data_path_Callback(hObject, eventdata, handles)
% hObject    handle to data_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of data_path as text
%        str2double(get(hObject,'String')) returns contents of data_path as a double


% --- Executes during object creation, after setting all properties.
function data_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function time_delay_Callback(hObject, eventdata, handles)
% hObject    handle to time_delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_delay as text
%        str2double(get(hObject,'String')) returns contents of time_delay as a double


% --- Executes during object creation, after setting all properties.
function time_delay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_delay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function time_trigger_Callback(hObject, eventdata, handles)
% hObject    handle to time_trigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_trigger as text
%        str2double(get(hObject,'String')) returns contents of time_trigger as a double


% --- Executes during object creation, after setting all properties.
function time_trigger_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_trigger (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in freq_modify.
function freq_modify_Callback(hObject, eventdata, handles)
% hObject    handle to freq_modify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of freq_modify



function fs_V_Callback(hObject, eventdata, handles)
% hObject    handle to fs_V (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fs_V as text
%        str2double(get(hObject,'String')) returns contents of fs_V as a double


% --- Executes during object creation, after setting all properties.
function fs_V_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs_V (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function fs_IQ_Callback(hObject, eventdata, handles)
% hObject    handle to fs_IQ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fs_IQ as text
%        str2double(get(hObject,'String')) returns contents of fs_IQ as a double


% --- Executes during object creation, after setting all properties.
function fs_IQ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs_IQ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function T_period_Callback(hObject, eventdata, handles)
% hObject    handle to T_period (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_period as text
%        str2double(get(hObject,'String')) returns contents of T_period as a double


% --- Executes during object creation, after setting all properties.
function T_period_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_period (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function T_scan_Callback(hObject, eventdata, handles)
% hObject    handle to T_scan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T_scan as text
%        str2double(get(hObject,'String')) returns contents of T_scan as a double


% --- Executes during object creation, after setting all properties.
function T_scan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T_scan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function filter_ICRF_Callback(hObject, eventdata, handles)
% hObject    handle to filter_ICRF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filter_ICRF as text
%        str2double(get(hObject,'String')) returns contents of filter_ICRF as a double


% --- Executes during object creation, after setting all properties.
function filter_ICRF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filter_ICRF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in inversion_method.
function inversion_method_Callback(hObject, eventdata, handles)
% hObject    handle to inversion_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns inversion_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from inversion_method


% --- Executes during object creation, after setting all properties.
function inversion_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inversion_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function freq_prob_Callback(hObject, eventdata, handles)
% hObject    handle to freq_prob (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of freq_prob as text
%        str2double(get(hObject,'String')) returns contents of freq_prob as a double


% --- Executes during object creation, after setting all properties.
function freq_prob_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freq_prob (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in judge_cutoff.
function judge_cutoff_Callback(hObject, eventdata, handles)
% hObject    handle to judge_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of judge_cutoff


% --- Executes on selection change in fitting_method.
function fitting_method_Callback(hObject, eventdata, handles)
% hObject    handle to fitting_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns fitting_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fitting_method


% --- Executes during object creation, after setting all properties.
function fitting_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fitting_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function rand_number_Callback(hObject, eventdata, handles)
% hObject    handle to rand_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rand_number as text
%        str2double(get(hObject,'String')) returns contents of rand_number as a double


% --- Executes during object creation, after setting all properties.
function rand_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rand_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in Run.
function Run_Callback(hObject, eventdata, handles)
% hObject    handle to Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% read the setting from the panel and save them
shot_1=get(handles.shot,'String');
IT_1=get(handles.IT,'String');
time_1=get(handles.time,'String');
data_path_1=get(handles.data_path,'String');
time_delay_1=get(handles.time_delay,'String'); 
time_trigger_1=get(handles.time_trigger,'String'); 
fs_V_1=get(handles.fs_V,'String'); 
fs_IQ_1=get(handles.fs_IQ,'String'); 
T_period_1=get(handles.T_period,'String'); 
T_scan_1=get(handles.T_scan,'String'); 
filter_freq_Q_1=get(handles.filter_freq_Q,'String'); 
filter_freq_V_1=get(handles.filter_freq_V,'String'); 
filter_freq_W_1=get(handles.filter_freq_W,'String'); 
filter_ICRF_1=get(handles.filter_ICRF,'String');
edge_profile_1=num2str(get(handles.edge_profile,'value')); 
edge_profile_function_1=num2str(get(handles.edge_profile_function,'value'));
judge_cutoff_1=num2str(get(handles.judge_cutoff,'value')); 
freq_prob_1=get(handles.freq_prob,'String'); 
rand_number_check_1=num2str(get(handles.rand_number_check,'value')); 
rand_number_1=get(handles.rand_number,'String');
judge_Rzero_1=num2str(get(handles.judge_Rzero,'value')); 
R_zero_1=get(handles.R_zero,'String');

save('refl_init.mat','shot_1','IT_1','time_1','data_path_1','time_delay_1', ... 
       'time_trigger_1','fs_V_1','fs_IQ_1','T_period_1','T_scan_1','filter_freq_Q_1','filter_freq_V_1', ...
       'filter_freq_W_1','filter_ICRF_1','edge_profile_1','edge_profile_function_1','judge_cutoff_1',...
   'freq_prob_1','rand_number_check_1','rand_number_1','judge_Rzero_1','R_zero_1');
   
shot=str2num(shot_1);       
IT=str2num(IT_1);
time=str2num(time_1);
data_path=data_path_1;
time_delay=str2num(time_delay_1)*10^-3; % ms to s
rand_number_check=str2num(rand_number_check_1);
rand_number=str2num(rand_number_1);  % The number for error analysis
judge_cutoff=str2num(judge_cutoff_1); % this is to judge whether calculate the cuoff positions using the ne profile
edge_profile=str2num(edge_profile_1);% to judge whether needs an assumed density profile?
edge_profile_function=str2num(edge_profile_function_1);
freq_prob=str2num(freq_prob_1); % probing frequency before X4
time_trigger=str2num(time_trigger_1)*10^-3; % ms to s
fs_V=str2num(fs_V_1)*10^6;  % MHZ
fs_IQ=str2num(fs_IQ_1)*10^6;  % MHz
T_period=str2num(T_period_1)*10^-6; % s to us
T_scan=str2num(T_scan_1)*10^-6; % s to us
filter_freq_Q=str2num(filter_freq_Q_1);
filter_freq_V=str2num(filter_freq_V_1);
filter_freq_W=str2num(filter_freq_W_1);
filter_ICRF=str2num(filter_ICRF_1)*10^6;
judge_Rzero=str2num(judge_Rzero_1);
R_zero=str2num(R_zero_1);
% close all

%% main calculation program

% mdsconnect('mds.ipp.ac.cn');
% mdsopen('east',shot);
% IT_raw=mdsvalue('\FOCS_IT');
% Itt=mdsvalue('dim_of(\FOCS_IT)');
% % [~,ind] = min(Itt-time);
% % IT = abs(IT_raw(ind));
% IT = abs(mean(IT_raw(time-0.5<Itt&Itt<time+0.5)));
% mdsdisconnect;
% save('FOCS_IT.mat','IT')

period=5; %the average periods 
t1=time-time_delay; % this time delay is due to the delay of other diagnostics.
f=freq_prob*10^9;

ICRF_on_V = 0;
ICRF_on_W = 0;

[fq,tau,tau_err]=refl_tau_analysis(shot,t1,IT,time_trigger,fs_V,fs_IQ,...
    T_period,T_scan,period,data_path,filter_freq_Q,filter_freq_V,...
    filter_freq_W,filter_ICRF,ICRF_on_V,ICRF_on_W);

fq_tau=[fq' tau'];
size(fq_tau)
%save(['matfile\' num2str(shot) '_fq_tau.mat'],'fq_tau')

[R1,ne1]=refl_density(IT,rand_number,fq,tau,tau_err,rand_number_check,...
    edge_profile,edge_profile_function,judge_Rzero,R_zero);
if rand_number_check==0
    
    mean_R=R1;
    mean_n=ne1;
   figure(100)
   plot(mean_R,mean_n,'bd')
   hold on
   X=[mean_R' mean_n'];
   save('X.mat','X')
   
elseif rand_number_check==1
R=reshape(R1,1,[]);
ne=reshape(ne1,1,[]);
save('R.mat','R1')
save('ne.mat','ne1')
figure(100)
plot(R,ne,'v')
mean_R=mean(R1,2);
mean_n=mean(ne1,2);
err_R=std(R1,0,2);
err_n=std(ne1,0,2);
hold on
errorbare('d',mean_R,mean_n,err_R,err_n,'-rv')
X=[R;ne];
save('X.mat','X')
end 

legend([num2str(shot),',t = ',num2str(time),' s'])
xlabel('R [m]')
ylabel('n_{e} [10^{19}m^{-3}]')
title(['EAST discharge #', num2str(shot), ', n_{e} profile from Refl.'])
%xlim([1.87 2.38])
ylim([0 inf])
hold on
% figure
% plot(mean_R,mean_n,'rs')
% cssss
[R_fit,ne_fit,nped,nwidth,ngrad]=refl_fitting(mean_R,mean_n);
figure(100)
% plot(R_fit,ne_fit,'color','black','linewidth',2)
n=[nped nwidth*100 ngrad];
save('n.mat','n')
% ssss
% figure
% plot(mean_R,mean_n,'s')

if judge_cutoff==1 %% calculate the cutoff positions of freq_prob
   1;
   f
    Bt=(IT*1.7)./(4086*R_fit); 
   for i=1:length(f)     
      [R_cut(i),ne_cut(i)]=cutoff(f(i),R_fit,Bt,ne_fit*10^19);
   end
    R_cut
    ne_cut
    figure(100)
    plot(R_cut,ne_cut/10^19,'o','color','r','markersize',10,'linewidth',3)
end

disp(['IT = ',num2str(IT)])
%dR=(R_cut(2:end)-R_cut(1))'
% guidata(hObject,handles)    

function [x1,y1,nped,nwidth,maxgrad]=refl_fitting(R,n)
    %% initial assumed value of MTANH fitting
    ped_pos=2.305;  %pedestal position    
    h=2; %height
    w=0.02; %width
    slope1=0.1; %slope of the core part  
    slope2=-0.07;% slope of edge part
    a0=[h/2 h/2 slope1 ped_pos w slope2];
%     index=find(R>2.2&R<=2.34);
    index=find(R>1);
    R1=R(index);
    n1=n(index);
    a=nlinfit(R1,n1,@mtanh_2,a0);
    x1=linspace(2,2.35,100);
    y1=mtanh_2(a,x1);
    nped=a(1)+a(2);
    nwidth=2*a(5);
    grad=abs(diff(y1)./diff(x1));
    maxgrad=max(grad);
    
       
function y=mtanh_1(a,x)
       %A=a(1);
       %B=a(2);
       h=a(1);
       alpha=a(2);
       x_sys=a(3);
       w=a(4); 
       beta=a(5);
       %gamma=a(6);
       z=(x_sys-x)./w;
       mtanh_fun=((1+alpha*z+0*z.^2).*exp(z)-(1+0*z).*exp(-z))./(exp(z)+exp(-z));
       y=h/2*(mtanh_fun+1);      
       
function y=mtanh_2(a,x)
       A=a(1);
       B=a(2);
       alpha=a(3);
       x_sys=a(4);
       w=a(5); 
       beta=a(6);
       %gamma=a(7);
       z=(x_sys-x)./w;
       mtanh_fun=((1+alpha*z).*exp(z)-(1+beta*z).*exp(-z))./(exp(z)+exp(-z));
       y=A*mtanh_fun+B;             
       
       
function [R_cut,ne_cut]=cutoff(f,R,Bt,ne)
         [R,index]=sort(R);
         ne=ne(index);
         Bt=Bt(index);
         f_R=rightcut(Bt,ne);
         plot(R,f_R);
         fdiff=(f-f_R);
         ind1=find(fdiff<0);
         ind2=find(fdiff>0);
         ind = max(ind1);
         if (~isempty(ind1)&&~isempty(ind2))
            % temp = isempty(ind)   
             R1=linspace(R(ind-1),R(ind+1),400);
             B1=spline(R,Bt,R1);
             n1=spline(R,ne,R1);
             u=X_group_V(f,B1,n1);
             ind1=find(abs(imag(u))==0);
             R2=R1(ind1);
             n2=n1(ind1);
             f_R1=rightcut(B1(ind1),n1(ind1));
             absdiff1=abs(f-f_R1);
             ind2=find(absdiff1==min(absdiff1));
             R_cut=R2(ind2);
             ne_cut=n2(ind2);
         else
             R_cut = -1;
             ne_cut = -1;
         end

function f_R=rightcut(Bt,ne)
         e=1.602*10^(-19);
         epsilon0=8.854*10^(-12);
         me=9.1096*10^(-31);
         c=3.0*10^8;
         Wpe=sqrt(e^2*ne./(epsilon0*me));
         Wce=e*Bt./me;
         W_R=sqrt(Wce.^2/4+Wpe.^2)+Wce/2;
         f_R=W_R/(2*pi);

function u=X_group_V(f,Bt,ne) 
         e=1.602*10^(-19);
         epsilon0=8.854*10^(-12);
         me=9.1096*10^(-31);
         c=3.0*10^8;
         Wpe=sqrt(e^2*ne./(epsilon0*me));
         Wce=e*Bt./me;
         X=(Wpe.^2)./((2*pi*f).^2);
         Y=Wce./(2*pi*f);
         Nx=sqrt(((1-X).^2-Y.^2)./(1-X-Y.^2)); 
         S2=(1-X-Y.^2).^2;
         S3=(1-X-Y.^2).^2+X.*(Y.^2);
         u=(c*Nx.*S2)./S3;
         
function N=X_refraction_index(f,Bt,ne)
        e=1.602*10^(-19);
         epsilon0=8.854*10^(-12);
         me=9.1096*10^(-31);
         c=3.0*10^8;
         Wpe=sqrt(e^2*ne./(epsilon0*me));
         Wce=e*Bt./me;
         X=(Wpe.^2)./((2*pi*f).^2);
         Y=Wce./(2*pi*f);
         N1=1-(X.*(1-X))./(1-X-Y.^2);
         index=find(N1<0);
         N1(index)=0;
         N=sqrt(N1);
         
function ne=find_ne(f,Bt)
         e=1.602*10^(-19);
         epsilon0=8.854*10^(-12);
         me=9.1096*10^(-31);
         c=3.0*10^8;
         Wce=e*Bt./me;
         WR=2*pi*f;
         X=(WR-Wce/2).^2-Wce.^2/4;
         ne=epsilon0*me*X./e^2;
         
function [freq1,fb1,freq2,fb2]=fb_plasma(freq,dfreqdt,R,Bt,ne)
    L1=length(freq);
    fb_plasma=zeros(1,L1);
    i1=zeros(1,L1);
    R_max=max(R);
    for j=1:L1
        f=freq(j);
        [R_cut,ne_cut]=cutoff(f,R,Bt,ne);
       if(R_cut ==-1 && ne_cut == -1);
         % no cutoff. wave reflect from the wall
         R1=linspace(R(1),R_max,5000);
         Bt1=spline(R,Bt,R1);
         ne1=spline(R,ne,R1);
         i1(j)=j;
       else
          % there exists cutoff in the plasma
        R1=linspace(R_cut,R_max,5000);
        Bt1=spline(R,Bt,R1);
        ne1=spline(R,ne,R1);
       end
        u_plasma=X_group_V(f,Bt1,ne1);
        fb(j)=2*dfreqdt*abs(trapz(R1,1./u_plasma));       
    end
    i1
    if ~isempty(find(i1==1))
        index1=find(i1==0);
        index2=find(i1>0);
        freq1=freq(index2);
        fb1=fb(index2);
        freq2=freq(index1);
        fb2=fb(index1);
    else
        freq1=freq;
        freq2=freq;
        fb1=fb;
        fb2=fb;
    end

function [R1,Ne]=refl_density(IT,rand_number,fq,tau,tau_err,rand_number_check,edge_profile,edge_profile_function,judge_Rzero,R_zero)
        R_max=2.5;
        R_inwall=1.356;
        R=linspace(R_inwall,R_max,1000);
        Bt=(IT*1.7)./(4086*R);  % Bt profile
        e=1.602*10^(-19);
        me=9.1096*10^(-31);
        c=3.0*10^8;  % light velocity
        if edge_profile==1
           %%using assumed edge profile if the edge reflectometry band is unavailable
          if edge_profile_function==1 %MTANH function
              ped_pos=2.28;  %pedestal position    
              h=2; %height
              w=0.03; %width
              slope=0.01; %slope of the core part  
              a0=[h/2 h/2 slope ped_pos w];
              R0=(R_max+R_inwall)/2; % major radius of device center
              R_plus =linspace(R0,R_max,500);
              R_minus =linspace(R_inwall,R0,500);
              ne_plus=MTANH(a0,R_plus)*10^19;   
              
          end          
           ne_minus = zeros(1,length(R_plus)-1);  
           for k=1:length(R_plus)-1
             ne_minus(k) = ne_plus(length(R_plus)+1-k);
           end
             ne = [ne_minus ne_plus];   % global density profile
             R=[R_minus(1:end-1) R_plus];
             Bt=(IT*1.7)./(4086*R); 
             freq=linspace(8,14,50)*4*10^9; % scan frequency range
             freq_scan=5*10^3;  %scan frequency 20 kHz
             dfreqdt=abs(freq(end)-freq(1))*freq_scan; %df/dt
             [freq1,fb1,freq2,fb2]=fb_plasma(freq,dfreqdt,R,Bt,ne);
             fb_temp=[fb1 fb2];
             tau_temp=fb_temp./dfreqdt*10^9; %ns
             c=3*10^8;
             tau_vac=2*(R_max-R_inwall)/c*10^9;%ns
             tau_assume=tau_temp-tau_vac;
             figure(1000)
             plot(freq/10^9,tau_assume,'linewidth',2)
             
             fq_1=[freq/10^9 fq];
             tau_1=[tau_assume tau];
             tau_err_1=[0.1*tau_assume tau_err];
             [fq_1,index]=sort(fq_1);
             tau_1=tau_1(index);
             tau_err_1=tau_err_1(index);
             
             R_zero=input('Please input the Zero density layer R_zero=');
             Bt_n0=(IT*1.7)./(4086*R_zero); 
             e=1.602*10^(-19);
             epsilon0=8.854*10^(-12);
             me=9.1096*10^(-31);
             fq_zero=e*Bt_n0/(2*pi*me)/10^9 %at zero density, fR=fce
             tau_zero=pchip(fq_1,tau_1,fq_zero);
             index=find(fq_1>fq_zero);
             fq_final=[fq_zero fq_1(index)];
             tau_final=[tau_zero tau_1(index)];
             tau_final_err=[tau_zero*0.1 tau_err_1(index)];
             figure
             errorbar(fq_final,tau_final,tau_final_err,'s')
             
        elseif judge_Rzero==1
            Bt_n0=(IT*1.7)./(4086*R_zero); 
            e=1.602*10^(-19);
            epsilon0=8.854*10^(-12);
            me=9.1096*10^(-31);
            fq_zero=e*Bt_n0/(2*pi*me)/10^9 %at zero density, fR=fce
            [fq_1,index]=sort(fq);
            tau_1=tau(index);
            tau_err_1=tau_err(index);
            tau_zero=pchip(fq_1,tau_1,fq_zero);
            index1=find(fq_1>fq_zero);
             fq_final=[fq_zero fq_1(index1)];
             tau_final=[tau_zero tau_1(index1)];
             tau_err_final=[tau_zero*0 tau_err_1(index1)];
             figure(1010)
             errorbar(fq_final,tau_final,tau_err_final,'-s')
        end
        
        fq_zero=fq_final(1); %frequency (GHz) for zero density
        Bt_n0=2*pi*fq_zero*10^9*me/e;  %at zero density, fR=fce
        R_n0=(IT*1.7)./(4086*Bt_n0)  %position of zero density
        display(['the zero density layer is at R=' num2str(R_n0) ,' m'])
        tau_vac=2*(R_n0-R_inwall)/c;
        tau1=tau_vac*10^9+tau_final;
        tau1_err=tau_err_final;
        fq1=fq_final;
        figure(1010)
        hold on
        errorbar(fq1,tau1,tau1_err*0.5,'-rv')
        
        
        R1=zeros(length(tau1),rand_number);
        Ne=zeros(length(tau1),rand_number);
        if rand_number_check==0
           tau_temp=tau1;
           phase=cumtrapz(2*pi*fq1*10^9,tau_temp/10^9);
           figure(1001)
           plot(fq1,phase,'-o')
           hold on
           [R1,ne]=refl_inversion_2(fq_zero*10^9,R_n0,0,fq1(2:end)*10^9,phase(2:end),IT);
           Ne=ne/10^19;

        elseif rand_number_check==1    
            for i=1:rand_number
                err_coef=rand(1,length(tau1_err))*2-1;
                tau_temp=tau1+1*err_coef.*tau1_err;
                %figure(1000)
                %plot(fq,tau_temp,'o','color',[i,0,i]/rand_number)
                phase=cumtrapz(2*pi*fq1*10^9,tau_temp/10^9);
                figure(1001)
                plot(fq1,phase)
                hold on
                %if inversion_method==1 %% simple inversion
                   [R,ne]=refl_inversion_2(fq_zero*10^9,R_n0,0,fq1(2:end)*10^9,phase(2:end),IT);
                %elseif inversion_method==2  %% slow inversion, this inversion is very slow.
                %   [R,ne]=refl_inversion_1(fq_zero*10^9,R_n0,0,fq(2:end)*10^9,phase(2:end),IT);
                %end

                R1(:,i)=R;
                Ne(:,i)=ne/10^19;
            end         
        end
        
 function [x3,n3]=refl_inversion_2(f0,x0,n0,freq,phase,IT)
        %simple inversion
        %f0 the lowest frequency, x0,n0 are the position and cutoff density
        %for f0. freq is a vector including all the other frequencies
        %larger than f0 and phase is the corresponding phase.
        x1=[];
        f1=[];
        c=3.0*10^8;
        for k=1:length(phase) 
            x=[x0 x1];
            f=[f0 f1];
            Bt1=(IT*1.7)./(4086*x);
            n1=find_ne(f,Bt1);
            N_temp=X_refraction_index(freq(k),Bt1,n1);
            if length(x)==1 
               phi_temp=0;
            else 
               phi_temp=-trapz(x,(4*pi*freq(k)/c)*N_temp);% note the sequence of x

            end
        dphi=phase(k)-phi_temp;
        x_temp=x(end)-(2*dphi*c/(4*pi*freq(k)))/N_temp(end);
        x1=[x1 x_temp];
        f1=[f1 freq(k)];
        end
      x2=[x0 x1];
      f2=[f0 freq];
      Bt2=(IT*1.7)./(4086*x2);
      n2=find_ne(f2,Bt2);
      x3=[x0 (x2(1:end-1)+x2(2:end))/2];
      n3=abs([n0 (n2(1:end-1)+n2(2:end))/2] );
        
  function [x2,n2]=refl_inversion_1(f0,x0,n0,freq,phase,IT)
        % slow inversion
        %f0 the lowest frequency, x0,n0 are the position and cutoff density
        %for f0. freq is a vector including all the other frequencies
        %larger than f0 and phase is the corresponding phase.
        x1=[];
        f1=[];
        c=3.0*10^8;
        for k=1:length(phase)
            x=[x0 x1];
            f=[f0 f1];
            Bt1=(IT*1.7)./(4086*x);
            %Bt1=spline(R,Bt,x);
            n1=find_ne(f,Bt1);
            freq1=freq(k);
            if length(x)==1 
               phase_temp1=0;
            else
               x_temp1=linspace(min(x),max(x),500);
               n_temp1=pchip(x,n1,x_temp1);
               Bt_temp1=(IT*1.7)./(4086*x_temp1);
               N_temp1=X_refraction_index(freq1,Bt_temp1,n_temp1);
               phase_temp1=abs(trapz(x_temp1,(4*pi*freq1/c)*N_temp1)); 
            end
            dphi=phase(k)-phase_temp1;
            x2=refl_phase(dphi,freq1,x(end),n1(end),R,Bt);    
        x1=[x1 x2];
        f1=[f1 freq1];
      end
      x2=[x0 x1];
      f2=[f0 freq];
      Bt2=(IT*1.7)./(4086*x2);
      n2=find_ne(f2,Bt2);      
        
function [fq_final,tau_final,tau_final_err]=refl_tau_analysis(shot,time,IT,time_trigger,fs_V,fs_IQ,T_period,T_scan,period,filepath,filter_freq_Q,filter_freq_V,filter_freq_W,filter_ICRF,ICRF_on_V,ICRF_on_W)
        % This function is to get the delay time from the Reflectometry
        % measurement 
        ts_V=1/fs_V;   %sampling frequency of voltage
        ts_IQ=1/fs_IQ; %sampling frequency of IQ signal 
        fftpoint=128;
        if shot<=52706 % Q: VO3260X_08_86259 V: VO3262P_01_86261  W:VO3262P_01_85757
            Q_f=[8.034 14.087];
            V_f=[11.82 19.149];
            W_f=[11.832 18.647];
            IQ_sign=[-1 1 1];
        elseif shot>52706&shot<55446 % Q: VO3260X_08_86259 V: VO3262P_01_86261 W: VO3262P_01_86261
            Q_f=[8.034 14.087];
            V_f=[11.82 19.149];
            W_f=[11.82 19.149];
            IQ_sign=[1 -1 -1];
        elseif shot>=55446&shot<57439  % Q: VO3260X_08_86259 V: VO3262P_01_86549 W: VO3262P_01_86549
            Q_f=[8.034 14.087];
            V_f=[11.6925 19.08];
            W_f=[11.6925 19.08];
            IQ_sign=[1 -1 -1];
        elseif shot>=57439&shot<57447  % Q: VO3260X_08_86259 V: VO3262P_01_86549 W: VO3262P_01_86549
            Q_f=[47 55]/4;
            V_f=[11.6925 19.08];
            W_f=[11.6925 19.08];
            IQ_sign=[1 -1 -1];   
            V_IQ_delay=0*10^-6;
        elseif shot>=57447&shot<57457
            Q_f=[47 52]/4;
            V_f=[50 60]/4;
            W_f=[11.6925 19.08];
            IQ_sign=[1 -1 -1];  
            V_IQ_delay=0*10^-6;
        elseif shot>=57457&shot<57460
             Q_f=[47 52]/4;
            V_f=[49 70]/4;
            W_f=[11.6925 19.08];
            IQ_sign=[1 -1 -1];  
            V_IQ_delay=0*10^-6;
        elseif shot>=57460&shot<57510
            Q_f=[8.034 14.087];
            V_f=[11.6925 19.08];
            W_f=[11.6925 19.08];
            IQ_sign=[1 -1 -1];
            V_IQ_delay=0*10^-6;
        elseif shot>=59056&shot<=77065
            Q_f=[8.1406 14.2073]+0.1; %upper band
            V_f=[11.7798 19.05]+0.1; %upper band
            W_f=[11.8313 19.1936]+0.1;%upper band
            IQ_sign=[1 1 1];
            V_IQ_delay=0*10^-6;
         elseif shot>=77065
            Q_f=[8.1406 14.2073]+0.1; %upper band
            V_f=[11.8427 19.3105]+0.1; %upper band
            W_f=[11.8313 19.1936]+0.1;%upper band
            IQ_sign=[1 1 1]; 
            V_IQ_delay=0.2*10^-6;
        end
        
        %% vacuum data
        [t_vac_IQ,Q_vac_data,V_vac_data,W_vac_data]=read_QVWdata_1(filepath,shot,time_trigger+0.004,time_trigger,period,T_period,T_scan,fs_V,fs_IQ,IQ_sign,V_IQ_delay);
        f_vac_Q=linspace(min(Q_f),max(Q_f),length(Q_vac_data))*4; % Q-band frequency (GHz)
        f_vac_V=linspace(min(V_f),max(V_f),length(V_vac_data))*4; % V-band frequency (GHz)
        f_vac_W=linspace(min(W_f),max(W_f),length(W_vac_data))*6; % W-band frequency (GHz)
        figure
        [fq_Q,tau_vac_Q,tau_err_vac_Q]=refl_fb_analysis(Q_vac_data,fftpoint,fs_IQ,filter_freq_Q,t_vac_IQ(1,:),f_vac_Q,filter_ICRF,[]);
        figure
        [fq_V,tau_vac_V,tau_err_vac_V]=refl_fb_analysis(V_vac_data,fftpoint,fs_IQ,filter_freq_V,t_vac_IQ(1,:),f_vac_V,filter_ICRF,[]);
        figure
        [fq_W,tau_vac_W,tau_err_vac_W]=refl_fb_analysis(W_vac_data,fftpoint,fs_IQ,filter_freq_W,t_vac_IQ(1,:),f_vac_W,filter_ICRF,[]);

        %% plasma data
        [t_IQ,Q_data,V_data,W_data]=read_QVWdata_1(filepath,shot,time,time_trigger,period,T_period,T_scan,fs_V,fs_IQ,IQ_sign,V_IQ_delay);
        f_Q=linspace(min(Q_f),max(Q_f),length(Q_data))*4; % GHz
        f_V=linspace(min(V_f),max(V_f),length(V_data))*4; % GHz
        f_W=linspace(min(W_f),max(W_f),length(W_data))*6; % GHz
        figure(55)
        [fq_Q,tau_Q,tau_err_Q]=refl_fb_analysis(Q_data,fftpoint,fs_IQ,[],t_IQ(1,:),f_Q,filter_ICRF,1);
        set(gcf,'Pointer','fullcross')
        x=ginput(2);
        x=sort(x);
        %x=[34 48.8];
        index=find(fq_Q>x(1)&fq_Q<x(2));
        fq_Q=fq_Q(index);
        tau1_Q=(tau_Q(index)-tau_vac_Q(index)); % ns
        tau1_err_Q=sqrt(tau_err_Q(index).^2+tau_err_vac_Q(index).^2);
        
        figure(56)
        [fq_V,tau_V,tau_err_V]=refl_fb_analysis(V_data,fftpoint,fs_IQ,[],t_IQ(1,:),f_V,filter_ICRF,1);
        set(gcf,'Pointer','fullcross')
        x=ginput(2);
        x=sort(x);
%         x=[49 72.2];
        index=find(fq_V>x(1)&fq_V<x(2))
                %%去除某些点
%         ind_2=find(fq_V>60&fq_V<60.5)
%         index = setdiff(index,ind_2)
%         ind_3=find(fq_V>61&fq_V<63.5)
%         index = setdiff(index,ind_3)
        fq_V=fq_V(index);
        
                
%         tau1_V=(tau_V(index)-tau_vac_V(index)); % ns
        
%         ICRF_on_V=0;
        if ICRF_on_V==1
            dfdt = (V_f(2)-V_f(1))*4*10^9/(40*10^-6);
            tau1_V = (tau_V(index)-tau_vac_V(index))+((34*10^6)./dfdt)*10^9; %ns    
        else 
            tau1_V=(tau_V(index)-tau_vac_V(index)); % ns
        end

        tau1_err_V=sqrt(tau_err_V(index).^2+tau_err_vac_V(index).^2);
        
        figure(57)
        [fq_W,tau_W,tau_err_W]=refl_fb_analysis(W_data,fftpoint,fs_IQ,[],t_IQ(1,:),f_W,filter_ICRF,1);
        set(gcf,'Pointer','fullcross')
        x=ginput(2);
        x=sort(x);
        %x=[1 2]
        index=find(fq_W>x(1)&fq_W<x(2));
        fq_W=fq_W(index);
        
        if ICRF_on_W==1
            dfdt = (W_f(2)-W_f(1))*4*10^9/(40*10^-6);
            tau1_W = (tau_W(index)-tau_vac_W(index))+((34*10^6)./dfdt)*10^9; %ns    
        else 
            tau1_W=(tau_W(index)-tau_vac_W(index)); % ns
        end
        tau1_err_W=sqrt(tau_err_W(index).^2+tau_err_vac_W(index).^2);

        figure(1000)
        errorbar(fq_Q,tau1_Q,tau1_err_Q,'o')
        hold on
        errorbar(fq_V,tau1_V,tau1_err_V,'rs')
        errorbar(fq_W,tau1_W,tau1_err_W,'v','color','black')
        %plot(fq,tau1)

       fq_final=[fq_Q fq_V fq_W];
       tau_final=[tau1_Q tau1_V tau1_W];
       tau_final_err=[tau1_err_Q tau1_err_V tau1_err_W];
       [fq_final,index1]=sort(fq_final);
       tau_final=tau_final(index1);
       
       %index=find(fq_final>=51.5&fq_final<=59);
       %tau_final(index)=[3.5598    3.1292   -0.1474   -3.1568   -4.6188   -5.1686   -5.5151   -5.7607   -5.7608];
       
       figure
       errorbar(fq_final,tau_final,tau_final_err,'-s')
        
function filt_data=fft_bandpass(data,band_pass,fs)
    %this function is to use FFT to filter the complex signal,
    %bandpass=filter_freq
    band_pass=sort(band_pass);
    L=length(data);
    y_data=fft(data);
    f=((1:L)*(fs/L)-fs/2)/10^6; % MHz
    y_data=fftshift(y_data);
    index=find(f<band_pass(1)|f>band_pass(2));
    y_data(index)=0;
    y_data=ifftshift(y_data);
    filt_data=ifft(y_data);

function [freq_m,tau,tau_err,filter_freq]=refl_fb_analysis(data,fftpoint,fs,filter_freq,t,f_band,filter_ICRF,judge_choose)
        s=size(data);
        %freq=linspace(freq_band(1),freq_band(2),s(2));
        g=1*fftpoint;
        m=0;
        %[A,B]=butter(2,filter_ICRF/(fs/2),'low'); % low pass to filter the ICRF frequency
        ind=[];
        while (m/2+1)*g<s(2)
                   k=(m*g/2+1):(m/2+1)*g;
                   ind=[ind;k];
                   x1=reshape(data(:,k)',1,[])';
                   %x1=filtfilt(A,B,x1);
                   freq1=f_band(k);
                   m=m+1;
                  [P1,f]=cpsd(x1-mean(x1),x1-mean(x1),[],[],fftpoint,fs);
                  P(:,m)=fftshift(P1);
                  freq_m(m)=mean(freq1);
%                   dfdt(m)=mean(abs(diff(freq1)./diff(t(k))))*10^12; 
                  dfdt(m)=(max(f_band)-min(f_band))*10^9./(40*10^-6);
        end
        
        f=f-f(fftpoint/2);
        f=f/10^6;
        
         imagesc(freq_m,f,log10(P));
         colormap('jet');
         set(gca,'YDir','normal')
        xlabel('frequency (GHz)','fontsize',14);
        ylabel('f_b(MHz)','fontsize',14)
        hold on
        
        
       if isempty(filter_freq)
           set(gcf,'Pointer','fullcross')
           [x100,y]=ginput(2);
           filter_freq=sort(y);
       else
           filter_freq=sort(filter_freq);
       end
       
       for i=1:m
          [fb(i),fb_err(i)]=find_fb(f,P(:,i),filter_freq);
       end
       
       %index=find(freq_m>=86&freq_m<=95);
       index=[];
       if isempty(index)|isempty(judge_choose)
       else
           for j=1:length(index)
               [fb(index(j)),fb_err(index(j))]=find_fb(f,P(:,index(j)),[7 14]);
           end
       end
       

       tau=((fb*10^6)./dfdt)*10^9; % ns
        tau_err=((fb_err*10^6)./dfdt)*10^9;
       errorbar(freq_m,fb,fb_err,'-o','color','k','markersize',6,...
           'linewidth',2)
       
       
 function [fb,fb_err]=find_fb(f,P,f_filter)
        index1=find(f>f_filter(1)&f<f_filter(2));
        f1=f(index1);
        P1=P(index1);
        index2=find(P1>=(max(P1)*exp(-1)));
        f2=f1(index2);
        P2=P1(index2);
        %fb=f2(find(P2==max(P2)));
       fb=sum(f2.*P2)/sum(P2);
        fb_err=sqrt(2*sum(((f2-fb).^2).*P2)/sum(P2));
       
function [t_IQ,Qdata_1,Vdata_1,Wdata_1]=read_QVWdata_1(filepath,shot,time,time_trigger,period,T_period,T_scan,fs_V,fs_IQ,IQ_sign,V_IQ_delay)
        %period: the number of period for average
        %T_period: the time of each period 
        %T_scan: the time for real frequency scan in one period T_period>T_scan
        %fs: the sampling frequency
        ts_V=1/fs_V; 
        ts_IQ=1/fs_IQ;
        period_length=str2num(num2str(T_period*fs_IQ))
        scan_length=str2num(num2str(T_scan*fs_IQ))
        period1=period*2;
       
        data_length_V=ceil(T_period*fs_V*period1)*2 %V/W voltage,read the Voltage data with a length of period1
        time_Voltage=(time+(0:(data_length_V/2-1))*ts_V)*1000;%ms
        
        t1=time-time_trigger% considering the delay between tigger time of Refl. sampling system  and EAST plasma
%         filename_Voltage=[filepath '\' num2str(shot) '\VCO.bin']; % the time series of QVW controlling voltage is the same, only V-band voltage is used.
%         fid_V= fopen(filename_Voltage);
%         %datapoint1=fix(fs_V*t1)*8*3;  % Q/V/W voltage, double precision
%         %datapoint1=str2num(num2str(datapoint1))
%         datapoint_V=(fs_V*t1*2)*2;% int16, two channels
%         datapoint_V=str2num(num2str(datapoint_V))
%         fseek(fid_V,datapoint_V,'cof');   %reposition the file position indicator
%         signal_Voltage= fread(fid_V,data_length_V,'int16');
%         %signal_Voltage=voltage_temp./max(voltage_temp)*20; 
%         signal_Voltage=reshape(signal_Voltage,2,[])/(0.65*10^4);
%         V_Voltage=signal_Voltage(1,:);
%         W_Voltage=signal_Voltage(2,:);
%         fclose(fid_V) 
%         figure
%         ax(1)=subplot(411)
%         plot(time_Voltage,V_Voltage,'r')
%         hold on
%         plot(time_Voltage,W_Voltage,'blue')
      %% read the I/Q data  
        filename_IQ=[filepath '\' num2str(shot) '\Q_V_W.bin'];  
        fid_IQ=fopen(filename_IQ);
        datapoint_IQ=(fs_IQ*t1*2)*8;% int16, 8 channels
        datapoint_IQ=str2num(num2str(datapoint_IQ))
        fseek(fid_IQ,datapoint_IQ,'cof');
        datalength_IQ=ceil(T_period*fs_IQ*period1)*8 %8 channels
        time_IQ=(time+(0:(datalength_IQ/8-1))*ts_IQ)*1000;%ms
        %datalength=str2num(num2str(datalength))
        data1=fread(fid_IQ,datalength_IQ,'int16');
        data1=reshape(data1,8,[]);
        Qdata=data1(1,:)+sqrt(-1)*IQ_sign(1)*data1(2,:);
        Vdata=data1(3,:)+sqrt(-1)*IQ_sign(2)*data1(4,:);
        Wdata=data1(5,:)+sqrt(-1)*IQ_sign(3)*data1(6,:);
        V_Voltage=data1(7,:)/(6*10^3);
        W_Voltage=data1(8,:)/(6*10^3);
        figure
         ax(1)=subplot(411)
         plot(time_Voltage,V_Voltage,'r')
         hold on
         plot(time_Voltage,W_Voltage,'blue')
        
        
        ax(2)=subplot(412)
        plot(time_IQ,real(Qdata))
        hold on
        plot(time_IQ,imag(Qdata),'r')
        ax(3)=subplot(413)
        plot(time_IQ,real(Vdata))
        hold on
        plot(time_IQ,imag(Vdata),'r')
        ax(4)=subplot(414)
        plot(time_IQ,real(Wdata))
        hold on
        plot(time_IQ,imag(Wdata),'r')
        linkaxes(ax,'x')
        index=find(V_Voltage==max(V_Voltage))
        index =index(1)    
        L=fix(index/period_length)
        ind=index-L*(period_length)
        subplot(411)
        plot(time_Voltage(ind),V_Voltage(ind),'ro')
        
        ind=ind+fix(V_IQ_delay*fs_IQ);       
        for i=1:period
            index1=ind+(i-1)*period_length+(period_length-scan_length);
            index=index1:(index1+scan_length-1);
            Qdata_1(i,:)=Qdata(index);
            Vdata_1(i,:)=Vdata(index);
            Wdata_1(i,:)=Wdata(index);
            t_IQ(i,:)=time_IQ(index);
        end
        
        

        fclose(fid_IQ)       
        
        figure
        ax1(1)=subplot(311);
        plot(t_IQ',real(Qdata_1)')
        hold on
        plot(t_IQ',imag(Qdata_1)','r')
        ax1(2)=subplot(312);
        plot(t_IQ',real(Vdata_1)')
        hold on
        plot(t_IQ',imag(Vdata_1)','r')
        ax1(3)=subplot(313);
        plot(t_IQ',real(Wdata_1)')
        hold on
        plot(t_IQ',imag(Wdata_1)','r')
        linkaxes(ax1,'x')
        
 function [t,data]=read_QVWdata_2(filename,time_trigger,time_select,period,T_period,T_scan,fs)
        %%read 2018 QVW data
        %period: the number of period for average
        %T_period: the time of each period 
        %T_scan: the time for real frequency scan in one period T_period>T_scan
        %fs: the sampling frequency
        ts=1/fs;
        fid=fopen(filename);
        time=time_select-period*T_period;  % s
        t1=time-time_trigger;
        datapoint=(fs*t1*2);% int16
        datapoint=str2num(num2str(datapoint));
        fseek(fid,datapoint,'cof');
        period1=period*3;
        datalength=fix(T_period*fs*period1);
        datalength=str2num(num2str(datalength));
        data1=fread(fid,datalength,'int16');
        t1=(time+(0:1:(datalength-1))*ts)*1000;%s to ms 
        dt=abs(t1-time_select*10^3);
        index=find(dt==min(dt));
        period_datapoint=str2num(num2str(T_period*fs));      
        scan_datapoint=str2num(num2str(T_scan*fs));     
        for i=1:period
            index1=index+(i-1)*period_datapoint;
            ind=index1:(index1+scan_datapoint-1);
            data(i,:)=data1(ind);
            t(i,:)=t1(ind);
        end
        fclose(fid)       
         
        
 function errorbare(sty,x,y,xbar,ybar,symbol)

% ERRORBARE Enhanced Errorbar Function.
%   ERRORBARE(STY,X,Y,Xbar,Ybar,symbol) 
%   It can draw errorbar along X/Y/Dual axis 
%   in normal,semilog,loglog coordinate system,
%   and adjust length of top line automatically,
%   can also control dotstyle and color in the same way with errorbar.
%
%   If the lower and upper error range of x/y is different, they should be
%   input as [lower,upper] if x/y is a column vector; 
%   for a row vector, they should be [lower;uper].
%
%   parameter STY include 12 types: 
%   v,h,d,vlogx,hlogx,dlogx,vlogy,hlogy,
%       dlogy,vlogd,hlogd,dlogd 
%   where
%   v stands for vertical errorbar，
%   h draws horizontal errorbar，
%   d means dual direction,
%   logx corresponding to semilogx，can use preffix v/h/d
%   logy corresponding to semilogy，can use preffix v/h/d
%   logd corresponding to loglog，can use preffix v/h/d

%   误差棒函数增强版
%   ERRORBARE(STY,X,Y,Xbar,Ybar,symbol) 
%   可在各个坐标系中沿X轴，Y轴方向，或者两轴方向绘制误差棒，
%   能够根据所选坐标类型调整端点线长。
%   增加对误差棒的线型控制，用法与原errorbar函数中相同
%
%   若上下限范围不同，X为列向量时应按照
%   [下限,上限] 的格式输入，若为行向量则为 [下限;上限]
%
%   STY 参数包括 v,h,d,vlogx,hlogx,dlogx,vlogy,hlogy,
%	dlogy,vlogd,hlogd,dlogd 共12种
%   v 表示误差棒垂直，
%   h 表示误差棒水平，
%   d (dual) 显示双轴误差，
%   logx 对应 semilogx，前缀 v,h,d 意义同上
%   logy 对应 semilogy，前缀 v,h,d 意义同上
%   logd 对应 loglog，前缀 v,h,d 意义同上

%   For example,
%	x = 1:10;
%	y = sin(x)+2;
%	e = std(y)*ones(size(x));
%	errorbare(x,y,e)	% use function "errorbar" directly
%	errorbare(x,y,e,'or')
%	errorbare('v',x,y,e)	% "e" is error of "y"
%	errorbare('v',x,y,[e;2*e])  % try different error limits
%	errorbare('hlogx',x,y,e)    % "e" is error of "x" here，
%	errorbare('d',x,y,e,e)
%	errorbare('d',x,y,e,e,'or')
%	errorbare('dlogd',x,y,e,e)
%
%   by Henry Sting  Email: henrysting@hotmail.com
%   $Revision: 1.2 $  $Date: 2010-2-8 $

lx=[];ux=[];ly=[];uy=[]; % 误差棒上下限
xl=[];xr=[];yl=[];yr=[]; % 端点短线左右限

if ~isstr(sty)
	if nargin == 3
		errorbar(sty,x,y)
		return
	elseif nargin == 4
		errorbar(sty,x,y,xbar)
		return
	elseif nargin > 4
		error('Please assign adopted coordinate system with symbol parameters.')
	end
elseif isstr(sty)
    if size(x)~=size(y)
		error('Coordinate array should be equal.')
	end
	if nargin == 4
		symbol='ob';
		if length(x)~=length(xbar)
			error('Format of Xbar is illegal.')
		end
		if sty(1) == 'v'
			ybar=xbar;xbar=[];
		elseif sty(1) == 'h'
			ybar=[];
		elseif sty(1) == 'd'
			error('Parameters are not enough.')
		else
			error('Symbol parameter is illegal.')
		end
	elseif nargin == 5 & ~isstr(ybar)
		symbol='ob';
		if length(x)~=length(xbar)
			error('Format of Xbar is illegal.')
		elseif length(y)~=length(ybar)
			error('Format of Ybar is illegal.')
		end
	elseif nargin == 5 & isstr(ybar)
		symbol=ybar;ybar=[];
		if length(x)~=length(xbar)
			error('Format of Xbar is illegal.')
		end
		if sty(1) == 'v'
			ybar=xbar;xbar=[];
		elseif sty(1) == 'h'
			ybar=[];
		elseif sty(1) == 'd'
			error('Parameters are not enough.')
		else
			error('Symbol parameter is illegal.')
		end
	elseif nargin == 6
		if length(x)~=length(xbar)
			error('Format of Xbar is illegal.')
		elseif length(y)~=length(ybar)
			error('Format of Ybar is illegal.')
		end
		if ~isstr(symbol)
			error('Symbol should be string')
		end
	end
end


[ls,col,mark,msg] = colstyle(symbol); if ~isempty(msg), error(msg); end
symbol = [ls mark col]; % Use marker only on data part
esymbol = ['-' col]; % Make sure bars are solid


% 转换为列距阵
[a,b]=size(x);
if a < b
	x=x';y=y';xbar=xbar';ybar=ybar';
	c=a;a=b;b=c;
end

%% 处理上下限不等
[xa,xb]=size(xbar);
if xb==1
	ux=xbar;lx=xbar; 
elseif xb==2
	lx=xbar(:,1);ux=xbar(:,2);
end

[ya,yb]=size(ybar);
if yb==1
	uy=ybar;ly=ybar; 
elseif yb==2
	ly=ybar(:,1);uy=ybar(:,2);
end

%% 描点
dx=(max(x(:))-min(x(:)))/100;
dy=(max(y(:))-min(y(:)))/100;
logn=10;
if length(sty) == 1
	xl = x-dx; xr = x+dx; yl = y-dy; yr = y+dy; % 定义端点短线长度
	plot(x,y,symbol);hold on
elseif length(sty) == 5 & sty(2:5) == 'logx' 
	dx=(log(max(x(:)))-log(min(x(:))))/100;
	xl = x/logn^dx;xr = x*logn^dx;yl = y-dy; yr = y+dy; 
	semilogx(x,y,symbol);hold on
elseif length(sty) == 5 & sty(2:5) == 'logy' 
	dy=(log(max(y(:)))-log(min(y(:))))/100;
	yl = y/logn^dy;yr = y*logn^dy;xl = x-dx; xr = x+dx; 
	semilogy(x,y,symbol);hold on
elseif length(sty) == 5 & sty(2:5) == 'logd' 
	dx=(log(max(x(:)))-log(min(x(:))))/100;
	dy=(log(max(y(:)))-log(min(y(:))))/100;
	xl = x/logn^dx;xr = x*logn^dx; yl = y/logn^dy;yr = y*logn^dy;
	loglog(x,y,symbol);hold on
end

%% 纵向
if sty(1) == 'v' | sty(1) == 'd'
vx = zeros(a*9,b);
vx(1:9:end,:) = x;
vx(2:9:end,:) = x;
vx(3:9:end,:) = NaN;
vx(4:9:end,:) = xl;
vx(5:9:end,:) = xr;
vx(6:9:end,:) = NaN;
vx(7:9:end,:) = xl;
vx(8:9:end,:) = xr;
vx(9:9:end,:) = NaN;

vy = zeros(a*9,b);
vy(1:9:end,:) = y-ly;
vy(2:9:end,:) = y+uy;
vy(3:9:end,:) = NaN;
vy(4:9:end,:) = y-ly;
vy(5:9:end,:) = y-ly;
vy(6:9:end,:) = NaN;
vy(7:9:end,:) = y+uy;
vy(8:9:end,:) = y+uy;
vy(9:9:end,:) = NaN;

plot(vx,vy,esymbol,'markersize',20)
end
%% 横向
if sty(1) == 'h' | sty(1) == 'd'
hx = zeros(a*9,b);
hx(1:9:end,:) = x-lx;
hx(2:9:end,:) = x+ux;
hx(3:9:end,:) = NaN;
hx(4:9:end,:) = x-lx;
hx(5:9:end,:) = x-lx;
hx(6:9:end,:) = NaN;
hx(7:9:end,:) = x+ux;
hx(8:9:end,:) = x+ux;
hx(9:9:end,:) = NaN;

hy = zeros(a*9,b);
hy(1:9:end,:) = y;
hy(2:9:end,:) = y;
hy(3:9:end,:) = NaN;
hy(4:9:end,:) = yl;
hy(5:9:end,:) = yr;
hy(6:9:end,:) = NaN;
hy(7:9:end,:) = yl;
hy(8:9:end,:) = yr;
hy(9:9:end,:) = NaN;

plot(hx,hy,esymbol,'markersize',20)
end      


% --- Executes on button press in rand_number_check.
function rand_number_check_Callback(hObject, eventdata, handles)
% hObject    handle to rand_number_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rand_number_check




% --- Executes on selection change in edge_profile_function.
function edge_profile_function_Callback(hObject, eventdata, handles)
% hObject    handle to edge_profile_function (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns edge_profile_function contents as cell array
%        contents{get(hObject,'Value')} returns selected item from edge_profile_function


% --- Executes during object creation, after setting all properties.
function edge_profile_function_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edge_profile_function (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end






% --- Executes on button press in edge_profile.
function edge_profile_Callback(hObject, eventdata, handles)
% hObject    handle to edge_profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edge_profile




% --- Executes on button press in judge_Rzero.
function judge_Rzero_Callback(hObject, eventdata, handles)
% hObject    handle to judge_Rzero (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of judge_Rzero



function R_zero_Callback(hObject, eventdata, handles)
% hObject    handle to R_zero (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R_zero as text
%        str2double(get(hObject,'String')) returns contents of R_zero as a double


% --- Executes during object creation, after setting all properties.
function R_zero_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R_zero (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



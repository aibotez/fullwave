function varargout = Refl_ray_tracing(varargin)
% REFL_RAY_TRACING M-file for Refl_ray_tracing.fig
%      REFL_RAY_TRACING, by itself, creates a new REFL_RAY_TRACING or raises the existing
%      singleton*.
%
%      H = REFL_RAY_TRACING returns the handle to a new REFL_RAY_TRACING or the handle to
%      the existing singleton*.
%
%      REFL_RAY_TRACING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REFL_RAY_TRACING.M with the given input arguments.
%
%      REFL_RAY_TRACING('Property','Value',...) creates a new REFL_RAY_TRACING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Refl_ray_tracing_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Refl_ray_tracing_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help Refl_ray_tracing

% Last Modified by GUIDE v2.5 27-Oct-2014 17:00:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Refl_ray_tracing_OpeningFcn, ...
                   'gui_OutputFcn',  @Refl_ray_tracing_OutputFcn, ...
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


% --- Executes just before Refl_ray_tracing is made visible.
function Refl_ray_tracing_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Refl_ray_tracing (see VARARGIN)

% Choose default command line output for Refl_ray_tracing
clc
fid=fopen('ray_tracing_init.mat');
if fid==-1
   shot_1='48929';
   IT_1='8000';
   time_1='3.2';
   frequency_1='60';
   w0_1='2';
   theta_1='5';
   RZ0_1='[4, -0.16]';
   w1_1='1';
   ray_number_1='10';
   choice_density_1='1';
   dtau_1='-0.000001';
   step_1='3000';
   save('ray_tracing_init.mat','shot_1','IT_1','time_1','frequency_1','w0_1','theta_1', ... 
       'RZ0_1','w1_1','ray_number_1','choice_density_1','dtau_1','step_1');
end
fclose all
load('ray_tracing_init.mat');
set(handles.shot,'string',shot_1);
set(handles.IT,'string',IT_1);
set(handles.time,'string',time_1);
set(handles.freq_prob,'string',frequency_1);
set(handles.beam_waist,'string',w0_1);
set(handles.cline_angle,'string',theta_1);
set(handles.center_coordinate,'string',RZ0_1);
set(handles.width_cal,'string',w1_1);
set(handles.ray_number,'string',ray_number_1);
set(handles.choice_density,'value',str2num(choice_density_1));
set(handles.dtau,'string',dtau_1);
set(handles.step,'string',step_1);

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Refl_ray_tracing wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Refl_ray_tracing_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function beam_waist_Callback(hObject, eventdata, handles)
% hObject    handle to beam_waist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of beam_waist as text
%        str2double(get(hObject,'String')) returns contents of beam_waist as a double


% --- Executes during object creation, after setting all properties.
function beam_waist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beam_waist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function cline_angle_Callback(hObject, eventdata, ~)
% hObject    handle to cline_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cline_angle as text
%        str2double(get(hObject,'String')) returns contents of cline_angle as a double


% --- Executes during object creation, after setting all properties.
function cline_angle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cline_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function center_coordinate_Callback(hObject, eventdata, handles)
% hObject    handle to center_coordinate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of center_coordinate as text
%        str2double(get(hObject,'String')) returns contents of center_coordinate as a double


% --- Executes during object creation, after setting all properties.
function center_coordinate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to center_coordinate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function width_cal_Callback(hObject, eventdata, handles)
% hObject    handle to width_cal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of width_cal as text
%        str2double(get(hObject,'String')) returns contents of width_cal as a double


% --- Executes during object creation, after setting all properties.
function width_cal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to width_cal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function ray_number_Callback(hObject, eventdata, handles)
% hObject    handle to ray_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ray_number as text
%        str2double(get(hObject,'String')) returns contents of ray_number as a double


% --- Executes during object creation, after setting all properties.
function ray_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ray_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
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


% --- Executes on selection change in choice_density.
function choice_density_Callback(hObject, eventdata, handles)
% hObject    handle to choice_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns choice_density contents as cell array
%        contents{get(hObject,'Value')} returns selected item from choice_density


% --- Executes during object creation, after setting all properties.
function choice_density_CreateFcn(hObject, eventdata, handles)
% hObject    handle to choice_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function dtau_Callback(hObject, eventdata, handles)
% hObject    handle to dtau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dtau as text
%        str2double(get(hObject,'String')) returns contents of dtau as a double


% --- Executes during object creation, after setting all properties.
function dtau_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dtau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function step_Callback(hObject, eventdata, handles)
% hObject    handle to step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of step as text
%        str2double(get(hObject,'String')) returns contents of step as a double


% --- Executes during object creation, after setting all properties.
function step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in ray_tracing_run.
function ray_tracing_run_Callback(hObject, eventdata, handles)
% hObject    handle to ray_tracing_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% read the setting from the panel and save them
shot_1=get(handles.shot,'String');
IT_1=get(handles.IT,'String');
time_1=get(handles.time,'String');
frequency_1=get(handles.freq_prob,'string');
w0_1=get(handles.beam_waist,'string');
theta_1=get(handles.cline_angle,'string');
RZ0_1=get(handles.center_coordinate,'string');
w1_1=get(handles.width_cal,'string');
ray_number_1=get(handles.ray_number,'string');
choice_density_1=num2str(get(handles.choice_density,'value'));
dtau_1=get(handles.dtau,'string');
step_1=get(handles.step,'string');

save('ray_tracing_init.mat','shot_1','IT_1','time_1','frequency_1','w0_1','theta_1', ... 
       'RZ0_1','w1_1','ray_number_1','choice_density_1','dtau_1','step_1');

shot=str2num(shot_1);
IT=str2num(IT_1);   
time=str2num(time_1);     
choice_density=str2num(choice_density_1);
%% 2D density distribution
if choice_density==1
   [FileName,PathName] = uigetfile('*.txt','Select the txt-file'); % 1D density profile from saved txt file
   ne=load([PathName,FileName]);
   s=size(ne);
   if s(1)>s(2)
       ne=ne';
   end
   R_n=ne(1,:);
   n=ne(2,:)*10^19;
   figure(1)
   plot(R_n,n/10^19)
   hold on
   xlabel('Major Radius R (m)','fontsize',14)
   ylabel('n_e (10^{19} m^{-3})','fontsize',14)
   title('density profile','fontsize',14)
   Z_n=zeros(size(R_n));
   [R_temp,Z,ne_temp]=EAST_2DnBt(R_n,Z_n,n,shot,time); %2D distribution of density, x is along R and y is along Z
   s = size(ne_temp);
   R_max = 5;
   dR = abs(R_temp(end)-R_temp(end-1));
   R1 = linspace((max(R_temp)+dR),R_max,100)';
   ne1 = zeros(length(R1),s(2));
   R = [R_temp;R1];
   ne = [ne_temp;ne1];   
   figure(2)
   contour(R,Z,ne'/10^19,50);
   colormap('jet');
   hold on
   xlabel('Major Radius R (m)','fontsize',14)
   ylabel('Vertical coordinate Z (m)','fontsize',14)
   
elseif choice_density==2
   [FileName,PathName] = uigetfile('*.mat','Select the mat-file'); % 2D density profile from saved mat file (R,Z,ne)
   load([FileName,PathName]);
   figure(2)
   contour(R,Z,ne'/10^19,50)
   hold on
   xlabel('Major Radius R (m)','fontsize',14)
   ylabel('Vertical coordinate Z (m)','fontsize',14)    
end

frequency=str2num(frequency_1)*10^9; %GHz
w0=str2num(w0_1)*10^(-2); % cm
w1=str2num(w1_1)*10^(-2); % cm 
theta=str2num(theta_1);
RZ0=str2num(RZ0_1);
ray_number=str2num(ray_number_1);
dtau=str2num(dtau_1);
step=str2num(step_1);
mode = 'O';
%% Initializing the reflectometry Gauss beam  
c=3*10^8; %light velocity
lamda=c/frequency;% wavelength in vacuum
k0=2*pi*frequency/c; % wave-number in vacuum
zR=pi*w0^2/lamda;
theta1=theta/180*pi; % the cline angle, conversion to radian
R0=RZ0(1) % the coordinate of the beam centre
Z0=RZ0(2)
d=linspace(0,w1,ray_number); % the distance of the intial points from the beam axis
x0=R0+d*sin(theta1);
y0=Z0+d*cos(theta1);

% % set the initial point from idea 1
k=k0+k0*((d.^2)./(2*zR^2)); 
kx0=-k*cos(theta1); %the initial kx
ky0=k*sin(theta1); %the initial ky
I0=-(d.^2)./(w0^2); % the initial value of I
P0=4*d.^2/(w0^4); %P=(ki)^2=(grad I)^2
beta_x0=8*(x0-R0)/w0^4;
beta_y0=8*(y0-Z0)/w0^4;
% % set the initial point from idea 2
% % x is z £¨R£©,y is x £¨Z£©,z is y £¨B£©
% yi = y0-Z0;
% zi = zeros(1,length(yi));
% R_cz = (zi.^2+zR^2)./zi.^2;
% R_cy = (yi.^2+zR^2)./yi.^2;
% w_z = sqrt(w0^2.*(1+zi.^2./zR.^2));
% w_y = sqrt(w0^2.*(1+yi.^2./zR.^2));
% kx0 = -sqrt(k0^2+4*yi.^2./w_y.^4)*cos(theta1);
% ky0 = sqrt(k0^2+4*yi.^2./w_y.^4)*sin(theta1);
% kz0 = k0.*zi./R_cz;
% 
% I0 = -(zi.^2./w_z.^2+yi.^2./w_y.^2);
% kiz0 = -2*zi./w_z.^2;
% kiy0 = -2*yi./w_y.^2;
% P0 = kiz0.^2+kiy0.^2;
% beta_x0 = 8.*zi./w_z.^4;
% beta_y0 = 8.*yi./w_y.^4;

tic
% for the central beam
[x1,y1]=ray_tracing(x0(1),y0(1),kx0(1),ky0(1),I0(1),P0(1),beta_x0(1),beta_y0(1),dtau,step,R,Z,ne,IT,frequency,2,mode);
% calculate the cutoff layer
R_judge=min(x1)+0.05;
% for the other ray
[x,y,kx,ky,P,beta_x,beta_y]=ray_tracing(x0,y0,kx0,ky0,I0,P0,beta_x0,beta_y0,dtau,step,R,Z,ne,IT,frequency,R_judge,mode);
toc
figure(2)
plot(x,y,'linewidth',1)    

for m=1:length(x0)
    x1=x(:,m);
    R_cut(m)=min(x1);
    index1=find(x1==min(x1));
    y_refl(m)=y(index1,m);
    x2=x(index1:end,m);
    y2=y(index1:end,m);
    dx=abs(x2-R0);
    index2=find(dx==min(dx));
    y_rece(m)=y2(index2);
end
figure(3)
plot(y_rece,y_refl,'-sb','linewidth',2)
hold on
ylabel('Z at cutoff (m)');
xlabel('Z at receiving plane (m)');
ind1=find(y_rece>-0.2&y_rece<1);
y_refl1=y_refl(ind1);
y_rece1=y_rece(ind1);
p=polyfit(y_rece1,y_refl1,1);
R_cutoff=mean(R_cut)
y_rece2 = -0.16:0.01:-0.12;
% y_rece2 = -0.125:0.01:-0.0644;
plot(y_rece2,polyval(p,y_rece2),'-rs','linewidth',2,'markersize',10)
y_refl_P1 = p(1)*(-0.12)+p(2);
y_refl_P2 = p(1)*(-0.16)+p(2);
line([-0.12,-0.12],[min(y_refl1),y_refl_P1],'linestyle','--','color','b','linewidth',2);
line([-0.16,-0.16],[min(y_refl1),y_refl_P2],'linestyle','--','color','b','linewidth',2);
dz = (y_refl_P2-y_refl_P1)*100
line([min(y_rece1),-0.12],[y_refl_P1,y_refl_P1],'linestyle','--','color','b','linewidth',2);
line([min(y_rece1),-0.16],[y_refl_P2,y_refl_P2],'linestyle','--','color','b','linewidth',2);
ylabel('Z at cutoff (m)');
xlabel('Z at receiving plane (m)');
% axis([-0.2,-0.08,-0.06,-0.02,]);
function [x,y,kx,ky,P,beta_x,beta_y]=ray_tracing(x0,y0,kx0,ky0,I0,P0,beta_x0,beta_y0,dtau,cycle_length,R,Z,ne,IT,freq,R_judge,mode)
% ray tracing calculation based on complex einkonal method for X-mode wave
% using RK4 numerical integration method
c=3*10^8; %light velocity
L=length(x0);
k0=2*pi*freq/c; %wave-number in vacuum
x=[];
y=[];
kx=[];
ky=[];
P=[];
beta_x=[];
beta_y=[];

i=1;
Dx=0.001;
Dy=0.001;

while i<cycle_length
    x=[x x0];
    y=[y y0];
    kx=[kx kx0];
    ky=[ky ky0];
    P=[P P0];
    beta_x=[beta_x beta_x0];
    beta_y=[beta_y beta_y0];
    x1=x0;
    y1=y0;
    kx1=kx0;
    ky1=ky0;
    P1=P0;
    beta_x1=beta_x0;
    beta_y1=beta_y0;
    outputs=XO_par(x1,y1,Dx,Dy,R,Z,ne,IT,freq,mode);
    [RK1_x,RK1_y,RK1_kx,RK1_ky]=RK4_XOmode(k0,kx1,ky1,P1,beta_x1,beta_y1,outputs,mode);
    
    x1=x0+RK1_x*dtau/2;
    y1=y0+RK1_y*dtau/2;
    kx1=kx0+RK1_kx*dtau/2;
    ky1=ky0+RK1_ky*dtau/2;
    outputs=XO_par(x1,y1,Dx,Dy,R,Z,ne,IT,freq,mode);
    [RK2_x,RK2_y,RK2_kx,RK2_ky]=RK4_XOmode(k0,kx1,ky1,P1,beta_x1,beta_y1,outputs,mode);
    
    x1=x0+RK2_x*dtau/2;
    y1=y0+RK2_y*dtau/2;
    kx1=kx0+RK2_kx*dtau/2;
    ky1=ky0+RK2_ky*dtau/2;
    outputs=XO_par(x1,y1,Dx,Dy,R,Z,ne,IT,freq,mode);
    [RK3_x,RK3_y,RK3_kx,RK3_ky]=RK4_XOmode(k0,kx1,ky1,P1,beta_x1,beta_y1,outputs,mode);
    
    x1=x0+RK3_x*dtau;
    y1=y0+RK3_y*dtau;
    kx1=kx0+RK3_kx*dtau;
    ky1=ky0+RK3_ky*dtau;
    outputs=XO_par(x1,y1,Dx,Dy,R,Z,ne,IT,freq,mode);
    [RK4_x,RK4_y,RK4_kx,RK4_ky]=RK4_XOmode(k0,kx1,ky1,P1,beta_x1,beta_y1,outputs,mode);
    
    RK_x=(RK1_x+2*RK2_x+2*RK3_x+RK4_x)/6;
    RK_y=(RK1_y+2*RK2_y+2*RK3_y+RK4_y)/6;
    RK_kx=(RK1_kx+2*RK2_kx+2*RK3_kx+RK4_kx)/6;
    RK_ky=(RK1_ky+2*RK2_ky+2*RK3_ky+RK4_ky)/6;
    x1_temp=x0+RK_x*dtau;
    y1_temp=y0+RK_y*dtau;
    kx1_temp=kx0+RK_kx*dtau;
    ky1_temp=ky0+RK_ky*dtau;
    
    if any(x1_temp<R_judge)
        P1_temp=zeros(1,L);
        beta_x0=zeros(1,L);
        beta_y0=zeros(1,L);
    else
        [k_ix_temp,k_iy_temp]=grad_I(x0,y0,I0,x1_temp,y1_temp,I0);
        P1_temp=k_ix_temp.^2+k_iy_temp.^2;
        [beta_x0,beta_y0]=grad_I(x0,y0,P0,x1_temp,y1_temp,P1_temp);
    end
    x0=x1_temp;
    y0=y1_temp;
    kx0=kx1_temp;
    ky0=ky1_temp;
    P0=P1_temp;
    i=i+1;
    if mod(i,200)==0
        disp(i)
    end
    if any([x0>=max(R),x0<=min(R),y0>=max(Z),y0<=min(Z)])
       break;   
    end
end

x=reshape([x x0],L,[])';
y=reshape([y y0],L,[])';
kx=reshape([kx kx0],L,[])';
ky=reshape([ky ky0],L,[])';
P=reshape([P P0],L,[])';
beta_x=reshape([beta_x beta_x0],L,[])';
beta_y=reshape([beta_y beta_y0],L,[])';


function [x,y,kx,ky,P,beta_x,beta_y]=ray_tracing_ver_old(x0,y0,kx0,ky0,I0,P0,beta_x0,beta_y0,dtau,cycle_length,R,Z,ne,IT,freq,R_judge)
%% ray tracing calculation based on complex einkonal method for X-mode wave
%% using RK4 numerical integration method
c=3*10^8; %light velocity
L=length(x0);
k0=2*pi*freq/c; %wave-number in vacuum
x=[];
y=[];
kx=[];
ky=[];
P=[];
beta_x=[];
beta_y=[];

i=1;
Dx=0.001;
Dy=0.001;
while i<cycle_length
    x=[x x0];
    y=[y y0];
    kx=[kx kx0];
    ky=[ky ky0];
    P=[P P0];
    beta_x=[beta_x beta_x0];
    beta_y=[beta_y beta_y0];
    x1=x0;
    y1=y0;
    kx1=kx0;
    ky1=ky0;
    P1=P0;
    beta_x1=beta_x0;
    beta_y1=beta_y0; 
    [f,df_dx,df_dy,dg_dx,dg_dy]=X_par(x1,y1,Dx,Dy,R,Z,ne,IT,freq);
    [RK1_x,RK1_y,RK1_kx,RK1_ky]=RK4_Omode(k0,kx1,ky1,P1,beta_x1,beta_y1,f,df_dx,df_dy,dg_dx,dg_dy);
     x1=x0+RK1_x*dtau/2;
     y1=y0+RK1_y*dtau/2;
     kx1=kx0+RK1_kx*dtau/2;
     ky1=ky0+RK1_ky*dtau/2;
     [f,df_dx,df_dy,dg_dx,dg_dy]=X_par(x1,y1,Dx,Dy,R,Z,ne,IT,freq);
     [RK2_x,RK2_y,RK2_kx,RK2_ky]=RK4_Omode(k0,kx1,ky1,P1,beta_x1,beta_y1,f,df_dx,df_dy,dg_dx,dg_dy);
     x1=x0+RK2_x*dtau/2;
     y1=y0+RK2_y*dtau/2;
     kx1=kx0+RK2_kx*dtau/2;
     ky1=ky0+RK2_ky*dtau/2;
     [f,df_dx,df_dy,dg_dx,dg_dy]=X_par(x1,y1,Dx,Dy,R,Z,ne,IT,freq);
     [RK3_x,RK3_y,RK3_kx,RK3_ky]=RK4_Omode(k0,kx1,ky1,P1,beta_x1,beta_y1,f,df_dx,df_dy,dg_dx,dg_dy);
     x1=x0+RK3_x*dtau;
     y1=y0+RK3_y*dtau;
     kx1=kx0+RK3_kx*dtau;
     ky1=ky0+RK3_ky*dtau;
     [f,df_dx,df_dy,dg_dx,dg_dy]=X_par(x1,y1,Dx,Dy,R,Z,ne,IT,freq);
     [RK4_x,RK4_y,RK4_kx,RK4_ky]=RK4_Omode(k0,kx1,ky1,P1,beta_x1,beta_y1,f,df_dx,df_dy,dg_dx,dg_dy);
     RK_x=(RK1_x+2*RK2_x+2*RK3_x+RK4_x)/6;
     RK_y=(RK1_y+2*RK2_y+2*RK3_y+RK4_y)/6;
     RK_kx=(RK1_kx+2*RK2_kx+2*RK3_kx+RK4_kx)/6;
     RK_ky=(RK1_ky+2*RK2_ky+2*RK3_ky+RK4_ky)/6;
     x1_temp=x0+RK_x*dtau;
     y1_temp=y0+RK_y*dtau;
     kx1_temp=kx0+RK_kx*dtau;
     ky1_temp=ky0+RK_ky*dtau;
     
     if x1_temp<R_judge
        P1_temp=zeros(1,L);
        beta_x0=zeros(1,L);
        beta_y0=zeros(1,L);
     else
     [k_ix_temp,k_iy_temp]=grad_I(x0,y0,I0,x1_temp,y1_temp,I0);
     P1_temp=k_ix_temp.^2+k_iy_temp.^2;
     [beta_x0,beta_y0]=grad_I(x0,y0,P0,x1_temp,y1_temp,P1_temp);
     end
     x0=x1_temp;
     y0=y1_temp;
     kx0=kx1_temp;
     ky0=ky1_temp;
     P0=P1_temp;
     i=i+1;
end
x=reshape([x x0],L,[])';
y=reshape([y y0],L,[])';
kx=reshape([kx kx0],L,[])';
ky=reshape([ky ky0],L,[])';
P=reshape([P P0],L,[])';
beta_x=reshape([beta_x beta_x0],L,[])';
beta_y=reshape([beta_y beta_y0],L,[])';





function [dI_dx,dI_dy]=grad_I(x0,y0,I0,x1,y1,I1)
  % the first point i=1 should be on the beam axis. Beam axis ray obey the
  % classical ray tracing, i.e. no need to consider the complex einkonal
   L=length(x0);
   dI_dx=zeros(1,L);
   dI_dy=zeros(1,L);
   for i=2:L;
        r1=[x0(i) y0(i)];
        r2=[x1(i) y1(i)];
        r3=[x1(i-1) y1(i-1)];
        [dI_dx(i),dI_dy(i)]=pd_3(r1,I0(i),r2,I1(i),r3,I1(i-1));
   end

    
function [df_dx,df_dy]=pd_3(x1,f1,x2,f2,x3,f3)
        %calculate the partial derivatives from three points
        df1=f2-f1;
        dx_1=x2(1)-x1(1); 
        dy_1=x2(2)-x1(2);
        df2=f3-f2;
        dx_2=x3(1)-x2(1);
        dy_2=x3(2)-x2(2);
        df_dx=(df2*dy_1-df1*dy_2)/(dx_2*dy_1-dx_1*dy_2);
        df_dy=(df2*dx_1-df1*dx_2)/(dy_2*dx_1-dy_1*dx_2);

        
function [RK_x,RK_y,RK_kx,RK_ky]=RK4_XOmode(k0,kx1,ky1,P1,beta_x1,beta_y1,outputs,mode)
switch mode
    
    case 'X'
        f = outputs(1,:);
        df_dx = outputs(2,:);
        df_dy = outputs(3,:);
        dg_dx = outputs(4,:);
        dg_dy = outputs(5,:);
        RK_x=-2*kx1.*f;
        RK_y=-2*ky1.*f;
        RK_kx=(-2*k0^2*f+kx1.^2+ky1.^2-P1).*df_dx+k0^2*dg_dx-beta_x1.*f;
        RK_ky=(-2*k0^2*f+kx1.^2+ky1.^2-P1).*df_dy+k0^2*dg_dy-beta_y1.*f;
        
    case 'O'
        dX_dx = outputs(1,:);
        dX_dy = outputs(2,:);
        RK_x=-2*kx1;
        RK_y=-2*ky1;
        RK_kx=k0^2*dX_dx-beta_x1;
        RK_ky=k0^2*dX_dy-beta_y1;
        
end

   
function [RK_x,RK_y,RK_kx,RK_ky]=RK4_Xmode(k0,kx1,ky1,P1,beta_x1,beta_y1,f,df_dx,df_dy,dg_dx,dg_dy) 
        RK_x=-2*kx1.*f;
        RK_y=-2*ky1.*f;
        RK_kx=(-2*k0^2*f+kx1.^2+ky1.^2-P1).*df_dx+k0^2*dg_dx-beta_x1.*f;
        RK_ky=(-2*k0^2*f+kx1.^2+ky1.^2-P1).*df_dy+k0^2*dg_dy-beta_y1.*f;
        
        
function [f,df_dx,df_dy,dg_dx,dg_dy]=X_par(x0,y0,Dx,Dy,R,Z,ne,IT,freq)
         e=1.602*10^(-19);
         epsilon0=8.854*10^(-12);
         me=9.1096*10^(-31);
         c=3.0*10^8;
         L=length(x0);
         x=[x0-Dx x0 x0 x0 x0+Dx];
         y=[y0 y0+Dy y0 y0-Dy y0];
         Bt=(IT*1.7)./(4086*x); % EAST toroidal magnetic field 
         index1=find(R>(min(x)-0.1)&R<(max(x)+0.1));
         index2=find(Z>(min(y)-0.1)&Z<(max(y)+0.1));
         n=interp2(R(index1),Z(index2),ne(index1,index2)',x,y,'spline'); 
         Wpe=sqrt(e^2*n./(epsilon0*me));
         Wce=e*Bt./me;
         X=(Wpe.^2)./((2*pi*freq).^2);
         Y=Wce./(2*pi*freq);
         f1=1-(X./(1-Y.^2));
         g1=((X.^2).*Y.^2)./(1-Y.^2).^2;
         f1=reshape(f1,L,[]);
         g1=reshape(g1,L,[]);
         f=f1(:,3)';
         df_dx=(((f1(:,5)-f1(:,3))./Dx+(f1(:,3)-f1(:,1))./Dx)*1/2)';
         dg_dx=(((g1(:,5)-g1(:,3))./Dx+(g1(:,3)-g1(:,1))./Dx)*1/2)';
         df_dy=(((f1(:,2)-f1(:,3))./Dy+(f1(:,3)-f1(:,4))./Dy)*1/2)';
         dg_dy=(((g1(:,2)-g1(:,3))./Dy+(g1(:,3)-g1(:,4))./Dy)*1/2)';
    
    function outputs=XO_par(x0,y0,Dx,Dy,R,Z,ne,IT,freq,mode)
        e=1.602*10^(-19);
        epsilon0=8.854*10^(-12);
        me=9.1096*10^(-31);
        c=3.0*10^8;
        L=length(x0);
        xx=[x0-Dx x0 x0 x0 x0+Dx];
        yy=[y0 y0+Dy y0 y0-Dy y0];
        Bt=(IT*1.7)./(4086*xx); % EAST toroidal magnetic field
        index1=find(R>(min(xx)-0.1)&R<(max(xx)+0.1));
        index2=find(Z>(min(yy)-0.1)&Z<(max(yy)+0.1));
        while length(index1)==1
            thres = 0.1*2;
            index1=find(R>(min(xx)-thres)&R<(max(xx)+thres));
        end
        while length(index2)==1
            thres = 0.1*2;
            index2=find(Z>(min(yy)-thres)&Z<(max(yy)+thres));
        end
        
        n=interp2(R(index1),Z(index2),ne(index1,index2)',xx,yy,'spline');
        Wpe=sqrt(e^2*n./(epsilon0*me));
        Wce=e*Bt./me;
        X=(Wpe.^2)./((2*pi*freq).^2);
        Y=Wce./(2*pi*freq);
        
        switch mode
            case 'X'
                f1=1-(X./(1-Y.^2));
                g1=((X.^2).*Y.^2)./(1-Y.^2).^2;
                f1=reshape(f1,L,[]);
                g1=reshape(g1,L,[]);
                f=f1(:,3)';
                df_dx=(((f1(:,5)-f1(:,3))./Dx+(f1(:,3)-f1(:,1))./Dx)*1/2)';
                dg_dx=(((g1(:,5)-g1(:,3))./Dx+(g1(:,3)-g1(:,1))./Dx)*1/2)';
                df_dy=(((f1(:,2)-f1(:,3))./Dy+(f1(:,3)-f1(:,4))./Dy)*1/2)';
                dg_dy=(((g1(:,2)-g1(:,3))./Dy+(g1(:,3)-g1(:,4))./Dy)*1/2)';
            case 'O'
                X = reshape(X,L,[]);
                dX_dx = (((X(:,5)-X(:,3))./Dx+(X(:,3)-X(:,1))./Dx)*1/2)';
                dX_dy = (((X(:,2)-X(:,3))./Dy+(X(:,3)-X(:,4))./Dy)*1/2)';
        end
        if mode == 'X'
            outputs=[f;df_dx;df_dy;dg_dx;dg_dy];
        elseif mode == 'O'
            outputs=[dX_dx;dX_dy];
        end
      
    function [f,g]=X_par_1(x0,y0,R,Z,ne,IT,freq)
          e=1.602*10^(-19);
          epsilon0=8.854*10^(-12);
          me=9.1096*10^(-31);
          c=3.0*10^8;
          Bt=(IT*1.7)./(4086*x0); % EAST toroidal magnetic field 
          n=interp2(R,Z,ne',x0,y0,'spline'); 
          Wpe=sqrt(e^2*n./(epsilon0*me));
          Wce=e*Bt./me;
          X=(Wpe.^2)./((2*pi*freq).^2);
          Y=Wce./(2*pi*freq);
          f=1-(X./(1-Y.^2));
          g=((X.^2).*Y.^2)./(1-Y.^2).^2;

function [r,z,n]=EAST_2DnBt(R,Z,ne,shot,time)
%input the density profile ne to acquire the 2D distribution of density
mdsconnect('mds.ipp.ac.cn'); %connect the mds service
mdsopen('efit_east',shot);
r=mdsvalue('_sig=\R');
z=mdsvalue('_sig=\Z');
psi=mdsvalue('_sig=\psirz');
psi_axis=mdsvalue('_sig=\ssimag');
psi_bdy=mdsvalue('_sig=\ssibry');
t=mdsvalue('dim_of(_sig)');
mdsdisconnect;
dt=abs(time-t);
index=find(dt==min(dt));
psi1=reshape(psi(:,:,index(1)),length(r),length(z));
psi1_axis=psi_axis(index(1));
psi1_bdy=psi_bdy(index(1));
psi1_norm=(psi1-psi1_axis)./(psi1_bdy-psi1_axis);
PSI=interp2(r,z,psi1_norm',R,Z,'spline');
s=size(psi1_norm);
psi2_norm=reshape(psi1_norm,1,[]);
n1=pchip(PSI,ne,psi2_norm);
index=find(psi2_norm>=max(PSI));

index = find(n1>8e19);
n1(index)=nan;
index = find(n1<0);

n1(index)=nan;
n=reshape(n1,s(1),s(2));






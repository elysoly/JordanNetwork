function varargout = GUI(varargin)
% GUI MATLAB code for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 07-Jun-2017 23:59:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
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


% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function hidden_neurons_Callback(hObject, eventdata, handles)
% hObject    handle to hidden_neurons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hidden_neurons as text
%        str2double(get(hObject,'String')) returns contents of hidden_neurons as a double


% --- Executes during object creation, after setting all properties.
function hidden_neurons_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hidden_neurons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function epsilon1_Callback(hObject, eventdata, handles)
% hObject    handle to epsilon1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epsilon1 as text
%        str2double(get(hObject,'String')) returns contents of epsilon1 as a double


% --- Executes during object creation, after setting all properties.
function epsilon1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epsilon1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% [data,timestamps] = xlsread('daily-minimum-temperatures-in-me.xlsx', 'Daily minimum temperatures in M', 'A:B');
% data=(data-min(data))/(max(data)-min(data));
data=handles.load_data.UserData;
l=length(data);
tr_data=data(1:round(0.7*l),:);
val_data=data(length(tr_data)+1:length(tr_data)+round(0.15*l),:);
l_prime=length(tr_data)+length(val_data)+1;
ts_data=data(l_prime:l);


epsilon1=str2num(handles.epsilon1.String);
epsilon2=str2num(handles.epsilon2.String);
hidden_neuron=str2num(handles.hidden_neurons.String);

[selected_a,selected_ar,selected_b,selected_bout,selected_z,selected_y]=jordan_net(tr_data,val_data,ts_data,epsilon1,epsilon2,hidden_neuron,handles)
% save net info

handles.run.UserData=[selected_a.'; selected_ar.';selected_b.';selected_bout.'; selected_z*ones(size(selected_a)).' ; selected_y.'; ];



function epsilon2_Callback(hObject, eventdata, handles)
% hObject    handle to epsilon2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epsilon2 as text
%        str2double(get(hObject,'String')) returns contents of epsilon2 as a double


% --- Executes during object creation, after setting all properties.
function epsilon2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epsilon2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)
cla(handles.axes1,'reset');
cla(handles.axes2,'reset');
cla(handles.axes3,'reset');
cla(handles.axes4,'reset');
cla(handles.axes5,'reset');
cla(handles.axes6,'reset');
cla(handles.axes7,'reset');
cla(handles.axes8,'reset');

% handles.hidden_layer_w_label.String='';
% handles.output_layer_w_label.String='';
% handles.nmi_label.String='0';
handles.test_board.String='Test measure: ';
clc;


% --- Executes on button press in save_net.
function save_net_Callback(hObject, eventdata, handles)

net=handles.run.UserData;
csvwrite('my_netwrok.dat',net);


% --- Executes on button press in test_net.
function test_net_Callback(hObject, eventdata, handles)
net=csvread('my_netwrok.dat');
selected_a=net(1,:).';
selected_ar=net(2,:).';
selected_b=net(3,:).';
selected_bout=net(4,:).';
selected_z=net(4,1);
selected_y=net(5,:).';

if(handles.mse.Value==1)
    measure_type='MSE';
end
if (handles.mae.Value==1)
        measure_type='MAE';
end
if (handles.rmae.Value==1)
    measure_type='RMAE';
end
if (handles.pi.Value==1)
    measure_type='PI';
end

data=handles.load_data.UserData;
l=length(data);
tr_data=data(1:round(0.7*l),:);
val_data=data(length(tr_data)+1:length(tr_data)+round(0.15*l),:);
l_prime=length(tr_data)+length(val_data)+1;
ts_data=data(l_prime:l);

axes(handles.axes4);
[ts_err,z_test]=test_jordan(measure_type,selected_a,selected_ar,selected_b,selected_bout,ts_data,selected_z,selected_y);
z_test(1)=rand(1);
plot(ts_data,'-r');
hold on 
plot(z_test,'-b');
legend('data','prediction');
ylabel('output and target');
title('Test data');
label=handles.test_board.String;
handles.test_board.String=strcat(label, num2str(ts_err));


% --- Executes on button press in load_data.
function load_data_Callback(hObject, eventdata, handles)
tmp=load('my_data.mat');
data=tmp.data;
handles.load_data.UserData=data;


% --- Executes on button press in batch_run.
function batch_run_Callback(hObject, eventdata, handles)

data=handles.load_data.UserData;
l=length(data);
tr_data=data(1:round(0.7*l),:);
val_data=data(length(tr_data)+1:length(tr_data)+round(0.15*l),:);
l_prime=length(tr_data)+length(val_data)+1;
ts_data=data(l_prime:l);


epsilon1=str2num(handles.epsilon1.String);
epsilon2=str2num(handles.epsilon2.String);
hidden_neuron=str2num(handles.hidden_neurons.String);

[selected_a,selected_ar,selected_b,selected_bout,selected_z,selected_y]=jordan_net_batch(tr_data,val_data,ts_data,epsilon1,epsilon2,hidden_neuron,handles)
% save net info

handles.run.UserData=[selected_a.'; selected_ar.';selected_b.';selected_bout.'; selected_z*ones(size(selected_a)).' ; selected_y.'; ];



function batch_size_Callback(hObject, eventdata, handles)
% hObject    handle to batch_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of batch_size as text
%        str2double(get(hObject,'String')) returns contents of batch_size as a double


% --- Executes during object creation, after setting all properties.
function batch_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to batch_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lr_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to lr_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lr_threshold as text
%        str2double(get(hObject,'String')) returns contents of lr_threshold as a double


% --- Executes during object creation, after setting all properties.
function lr_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lr_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pause_time_Callback(hObject, eventdata, handles)
% hObject    handle to pause_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pause_time as text
%        str2double(get(hObject,'String')) returns contents of pause_time as a double


% --- Executes during object creation, after setting all properties.
function pause_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pause_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

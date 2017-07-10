function varargout = Bayesian(varargin)
% BAYESIAN MATLAB code for Bayesian.fig
%      BAYESIAN, by itself, creates a new BAYESIAN or raises the existing
%      singleton*.
%
%      H = BAYESIAN returns the handle to a new BAYESIAN or the handle to
%      the existing singleton*.
%
%      BAYESIAN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BAYESIAN.M with the given input arguments.
%
%      BAYESIAN('Property','Value',...) creates a new BAYESIAN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Bayesian_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Bayesian_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Bayesian

% Last Modified by GUIDE v2.5 10-Jul-2017 14:44:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Bayesian_OpeningFcn, ...
                   'gui_OutputFcn',  @Bayesian_OutputFcn, ...
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

function Bayesian_OpeningFcn(hObject, eventdata, handles, varargin)
movegui('center') 
axes(handles.logo)
matlabImage = imread('physicslogo.png');
image(matlabImage)
axis off
axis image
handles.output = hObject;
guidata(hObject, handles);

function varargout = Bayesian_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
%---------------------------------Unused Callbacks-------------------------
function interval_list_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function xlim_Callback(hObject, eventdata, handles)
function xlim_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function window_size_Callback(hObject, eventdata, handles)
function window_size_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function overlap_Callback(hObject, eventdata, handles)
function overlap_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function prop_const_Callback(hObject, eventdata, handles)
function prop_const_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function order_Callback(hObject, eventdata, handles)
function order_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function length_Callback(hObject, eventdata, handles)
function length_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function intervals_Callback(hObject, eventdata, handles)
function intervals_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function status_Callback(hObject, eventdata, handles)
function status_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function freq_1_Callback(hObject, eventdata, handles)
function freq_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function freq_2_Callback(hObject, eventdata, handles)
function freq_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function signal_list_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function interval_list_ButtonDownFcn(hObject, eventdata, handles)
%---------------------------------Unused Callbacks end---------------------
function refresh_limits_Callback(hObject, eventdata, handles)
x = get(handles.time_series_1,'xlim');
t = x(2) - x(1);
x = strcat(num2str(x(1)),' , ',num2str(x(2)));    

set(handles.xlim,'String',x);
set(handles.length,'String',t);

function add_interval_Callback(hObject, eventdata, handles)
f1 = str2double(get(handles.freq_1,'String'));
f2 = str2double(get(handles.freq_2,'String'));
ws = str2double(get(handles.window_size,'String'));
if isnan(ws)
    ws = 1/min(f1,f2);
end
fl = sprintf('%.3f,%.3f | %.2f',min(f1,f2),max(f1,f2),ws);
list = get(handles.interval_list,'String');
list{end+1,1} = fl;
set(handles.interval_list,'String',list);      
    
function interval_list_KeyPressFcn(hObject, eventdata, handles)
switch eventdata.Key
    case 'delete'
        interval_selected = get(handles.interval_list,'Value');
        if min(interval_selected)>1
            set(handles.interval_list,'Value',min(interval_selected)-1);
        else
            set(handles.interval_list,'Value',1);
        end
        list = get(handles.interval_list,'String');
        list(interval_selected,:) = [];
        set(handles.interval_list,'String',list);
        handles.bands(:,interval_selected) = [];
        guidata(hObject,handles);
        interval_list_Callback(hObject, eventdata, handles)
        guidata(hObject,handles);
        drawnow;
end   

%---------------------------Reading the data
function file_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function csv_read_Callback(hObject, eventdata, handles)
%Read csv file
    linkaxes([handles.time_series_1 handles.time_series_2],'x');
    set(handles.status,'String','Importing Signal...');
    handles.fs = str2double(cell2mat(inputdlg('Enter the Sampling Frequency')));
    fs = handles.fs;
    if isnan(fs)
      errordlg('Sampling frequency must be specified before importing','Parameter Error');
      return;
    end
    
    sig = read_from_csv();
    
    % Construct a questdlg with three options
    choice = questdlg('Select Orientation of Data set?', ...
        'Data Import','Column wise','Row wise','default');
    switch choice
        case 'Column wise'
            sig = sig';
            %msgbox(sprintf('Number of data sets %d',size(sig,1)),'Import Complete');
        case 'Row wise'
            %msgbox(sprintf('Number of data sets %d',size(sig,1)),'Import Complete');
        case 'default'
            errordlg('Data set orientation must be specified')
            return;
    end
    list = cell(size(sig,1)/2+1,1);
    list{1,1} = 'Signal Pair 1';
    
    for i = 2:size(sig,1)/2
        list{i,1} = sprintf('Signal Pair %d',i);
    end
    set(handles.signal_list,'String',list);

    handles.sig = sig;   
    time = 1:size(sig,2);
    time = time./fs;
    handles.time_axis = time;
    
    plot(handles.time_series_1,time,sig(1,:));%Plotting the time_series part after calculation of appropriate limits
    xlim(handles.time_series_1,[0,size(sig,2)./fs]);
    plot(handles.time_series_2,time,sig(1+size(sig,1)/2,:));%Plotting the time_series part after calculation of appropriate limits
    xlim(handles.time_series_2,[0,size(sig,2)./fs]);
    
    refresh_limits_Callback(hObject, eventdata, handles);%updates the values in the box
    guidata(hObject,handles);  
    xlabel(handles.time_series_2,'Time (s)');
    set(handles.status,'String','Select Data And Continue With Wavelet Transform');

% --------------------------------------------------------------------
function mat_read_Callback(hObject, eventdata, handles)
%Read mat file    
    linkaxes([handles.time_series_1 handles.time_series_2],'x');
    set(handles.status,'String','Importing Signal...');
    handles.fs = str2double(cell2mat(inputdlg('Enter the Sampling Frequency')));
    fs = handles.fs;
    if isnan(fs)
      errordlg('Sampling frequency must be specified before importing','Parameter Error');
      return;
    end
    
    sig = read_from_mat();
    sig = struct2cell(sig);
    sig = cell2mat(sig);
    
    % Construct a questdlg with three options
    choice = questdlg('Select Orientation of Data set?', ...
        'Data Import','Column wise','Row wise','default');
    switch choice
        case 'Column wise'
            sig = sig';
            %msgbox(sprintf('Number of data sets %d',size(sig,1)),'Import Complete');
        case 'Row wise'
            %msgbox(sprintf('Number of data sets %d',size(sig,1)),'Import Complete');
        case 'default'
            errordlg('Data set orientation must be specified')
            return;
    end
    list = cell(size(sig,1)/2+1,1);
    list{1,1} = 'Signal Pair 1';
    for i = 2:size(sig,1)/2
        list{i,1} = sprintf('Signal Pair %d',i);
    end
    set(handles.signal_list,'String',list);

    handles.sig = sig;   
    time = 1:size(sig,2);
    time = time./fs;
    handles.time_axis = time;
    
    plot(handles.time_series_1,time,sig(1,:));%Plotting the time_series part after calculation of appropriate limits
    xlim(handles.time_series_1,[0,size(sig,2)./fs]);
    plot(handles.time_series_2,time,sig(1+size(sig,1)/2,:));%Plotting the time_series part after calculation of appropriate limits
    xlim(handles.time_series_2,[0,size(sig,2)./fs]);
    
    refresh_limits_Callback(hObject, eventdata, handles);%updates the values in the box
    guidata(hObject,handles);  
    xlabel(handles.time_series_2,'Time (s)');
    set(handles.status,'String','Select Data And Continue With Wavelet Transform');



function calculate_Callback(hObject, eventdata, handles)
xl = get(handles.xlim,'String');
xl = csv_to_mvar(xl);
xl = xl.*handles.fs;
xl(2) = min(xl(2),size(handles.sig,2));
xl(1) = max(xl(1),1);
handles.sig_cut = handles.sig(:,xl(1):xl(2));

band_list = get(handles.interval_list,'String'); 
warning off;
tic
for i = 1:size(handles.sig,1)
    for j =1:size(band_list,1)
        fl = csv_to_mvar(band_list{j,1}(1:strfind(band_list{j,1},'|')-1));
        [handles.bands{i,j},~] = loop_butter(handles.sig_cut(i,:),fl,handles.fs);
    end
end
toc
guidata(hObject,handles);
warning on;
linkaxes([handles.phi1_axes,handles.phi2_axes,handles.coupling_strength_axis,...
    handles.time_series_1,handles.time_series_2],'x');
signal_list_Callback(hObject, eventdata, handles)

function signal_list_Callback(hObject, eventdata, handles)
%Selecting the signal pair   
signal_selected = get(handles.signal_list, 'Value');        
plot(handles.time_series_1, handles.time_axis, handles.sig(signal_selected,:));%Plotting the time_series part after calculation of appropriate limits
xl = csv_to_mvar(get(handles.xlim, 'String'));
xlim(handles.time_series_1, xl);
plot(handles.time_series_2, handles.time_axis, handles.sig(signal_selected+size(handles.sig,1)/2,:));%Plotting the time_series part after calculation of appropriate limits
xlim(handles.time_series_2, xl);        
xlabel(handles.time_series_2, 'Time (s)');
refresh_limits_Callback(hObject, eventdata, handles);%updates the values in the box 
set(handles.status, 'String', 'Select Data And Continue With Wavelet Transform');
if isfield(handles,'TPC')
    xyplot_Callback(hObject, eventdata, handles);
end
interval_list_Callback(hObject, eventdata, handles)
intervals_Callback(hObject, eventdata, handles)  
    
function interval_list_Callback(hObject, eventdata, handles)
cla(handles.coupling_strength_axis,'reset');
interval_selected = get(handles.interval_list,'Value');
signal_selected = get(handles.signal_list,'Value');
ovr = str2double(get(handles.overlap,'String'));
pr = str2double(get(handles.prop_const,'String'));
bn =  str2double(get(handles.order,'String'));
win = str2double(get(handles.window_size,'String'));%FIX THIS
tic
for i = 1:size(interval_selected,2)
    phi1 = angle(hilbert(handles.bands{signal_selected,interval_selected(i)}));          
    plot(handles.phi1_axes, handles.time_axis, phi1);
    phi2 = angle(hilbert(handles.bands{signal_selected+size(handles.sig,1)/2,interval_selected(i)}));         
    plot(handles.phi2_axes, handles.time_axis, phi2);
    [tm,cc,~] = bayes_main(phi1,phi2,win,1/handles.fs,ovr,pr,0,bn);
    for j = 1:size(cc,1)
        [cpl1(j),cpl2(j),~] = dirc(cc(j,:),bn);
    end
    hold(handles.coupling_strength_axis,'on');
    plot(handles.coupling_strength_axis,tm,cpl1); 
    plot(handles.coupling_strength_axis,tm,cpl2); 
    legend(handles.coupling_strength_axis,'cp1','cp2');
end 
toc




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

% Last Modified by GUIDE v2.5 07-Jul-2017 18:34:31

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


function interval_list_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function xlim_Callback(hObject, eventdata, handles)

function xlim_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function length_Callback(hObject, eventdata, handles)

function length_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function refresh_limits_Callback(hObject, eventdata, handles)
x = get(handles.time_series_1,'xlim');
t = x(2) - x(1);
x = strcat(num2str(x(1)),' , ',num2str(x(2)));    

set(handles.xlim,'String',x);
set(handles.length,'String',t);

function add_interval_Callback(hObject, eventdata, handles)
f1 = str2double(get(handles.freq_1,'String'));
f2 = str2double(get(handles.freq_2,'String'));
fl = sprintf('%f,%f',f1,f2);
list = get(handles.interval_list,'String');
list{end+1,1} = fl;
set(handles.interval_list,'String',list);

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

function signal_list_1_Callback(hObject, eventdata, handles)
signal_selected = get(handles.signal_list_1, 'Value');
plot(handles.time_series_1, handles.time_axis, handles.sig1(signal_selected,:));%Plotting the time_series part after calculation of appropriate limits
xl = csv_to_mvar(get(handles.xlim, 'String'));
xlim(handles.time_series_1, xl);
xlabel(handles.time_series_1, 'Time (s)');
refresh_limits_Callback(hObject, eventdata, handles);%updates the values in the box
xlabel(handles.time_series_1, 'Time (s)');
set(handles.status, 'String', 'Select Data And Continue With Wavelet Transform');
function signal_list_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function signal_list_2_Callback(hObject, eventdata, handles)
signal_selected = get(handles.signal_list_2, 'Value');
plot(handles.time_series_2, handles.time_axis, handles.sig1(signal_selected,:));%Plotting the time_series part after calculation of appropriate limits
xl = csv_to_mvar(get(handles.xlim, 'String'));
xlim(handles.time_series_2, xl);
xlabel(handles.time_series_2, 'Time (s)');
refresh_limits_Callback(hObject, eventdata, handles);%updates the values in the box
xlabel(handles.time_series_2, 'Time (s)');
set(handles.status, 'String', 'Select Data And Continue With Wavelet Transform');
function signal_list_2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function interval_list_Callback(hObject, eventdata, handles)
function interval_list_ButtonDownFcn(hObject, eventdata, handles)
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
        drawnow;
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

%---------------------------Reading the data
function file_Callback(hObject, eventdata, handles)
function sig1_Callback(hObject, eventdata, handles)
function sig2_Callback(hObject, eventdata, handles)

function csv_read1_Callback(hObject, eventdata, handles)
set(handles.status,'String','Importing Signal...');
    handles.fs = str2double(cell2mat(inputdlg('Enter the Sampling Frequency')));    
    
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
    
    list = cell(size(sig,1)+1,1);
    list{1,1} = 'Signal 1';
    for i = 2:size(sig,1)
        list{i,1} = sprintf('Signal %d',i);
    end
    set(handles.signal_list_1,'String',list);

    time = 1:size(sig,2);
    time = time./handles.fs;
    handles.time_axis = time;
    plot(handles.time_series_1,time,sig(1,:));%Plotting the time_series part after calculation of appropriate limits
    xlim(handles.time_series_1,[0,size(sig,2)./handles.fs]);
    
    handles.sig1 = sig;
    
    refresh_limits_Callback(hObject, eventdata, handles);%updates the values in the box
    guidata(hObject,handles);  
    set(handles.status,'String','Select Data And Continue With Wavelet Transform');
    
function mat_read1_Callback(hObject, eventdata, handles)
set(handles.status,'String','Importing Signal...');
    handles.fs = str2double(cell2mat(inputdlg('Enter the Sampling Frequency')));   
    
    sig = cell2mat(struct2cell(read_from_mat()));
    
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
    
    list = cell(size(sig,1)+1,1);
    list{1,1} = 'Signal 1';
    for i = 2:size(sig,1)
        list{i,1} = sprintf('Signal %d',i);
    end
    set(handles.signal_list_1,'String',list);

    time = 1:size(sig,2);
    time = time./handles.fs;
    handles.time_axis = time;
    plot(handles.time_series_1,time,sig(1,:));%Plotting the time_series part after calculation of appropriate limits
    xlim(handles.time_series_1,[0,size(sig,2)./handles.fs]);
    
    handles.sig1 = sig;
    
    refresh_limits_Callback(hObject, eventdata, handles);%updates the values in the box
    guidata(hObject,handles);  
    set(handles.status,'String','Select Data And Continue With Wavelet Transform');

function csv_read2_Callback(hObject, eventdata, handles)
set(handles.status,'String','Importing Signal...');
    handles.fs = str2double(cell2mat(inputdlg('Enter the Sampling Frequency'))); 

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
    
    list = cell(size(sig,1)+1,1);
    list{1,1} = 'Signal 1';
    for i = 2:size(sig,1)
        list{i,1} = sprintf('Signal %d',i);
    end
    set(handles.signal_list_2,'String',list);
      
    time = 1:size(sig,2);
    time = time./handles.fs;
    handles.time_axis = time;
    plot(handles.time_series_2,time,sig(1,:));%Plotting the time_series part after calculation of appropriate limits
    xlim(handles.time_series_2,[0,size(sig,2)./handles.fs]);
    xlabel(handles.time_series,'Time (s)');
    handles.sig2 = sig;
    refresh_limits_Callback(hObject, eventdata, handles);%updates the values in the box
    guidata(hObject,handles);  
    set(handles.status,'String','Select Data And Continue With Wavelet Transform');
    
function mat_read2_Callback(hObject, eventdata, handles)
set(handles.status,'String','Importing Signal...');
    handles.fs = str2double(cell2mat(inputdlg('Enter the Sampling Frequency')));

    sig = cell2mat(struct2cell(read_from_mat()));
    
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
    
    list = cell(size(sig,1)+1,1);
    list{1,1} = 'Signal 1';
    for i = 2:size(sig,1)
        list{i,1} = sprintf('Signal %d',i);
    end
    set(handles.signal_list_2,'String',list);
      
    time = 1:size(sig,2);
    time = time./handles.fs;
    handles.time_axis = time;
    plot(handles.time_series_2,time,sig(1,:));%Plotting the time_series part after calculation of appropriate limits
    xlim(handles.time_series_2,[0,size(sig,2)./handles.fs]);
    xlabel(handles.time_series_2,'Time (s)');
    handles.sig2 = sig;
    refresh_limits_Callback(hObject, eventdata, handles);%updates the values in the box
    guidata(hObject,handles);  
    set(handles.status,'String','Select Data And Continue With Wavelet Transform');


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

function calculate_Callback(hObject, eventdata, handles)

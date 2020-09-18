function varargout = phoebe(varargin)
%DESCRIPTION

%PHOEBE MATLAB code for phoebe.fig
%      PHOEBE, by itself, creates a new PHOEBE or raises the existing
%      singleton*.
%
%      H = PHOEBE returns the handle to a new PHOEBE or the handle to
%      the existing singleton*.
%
%      PHOEBE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PHOEBE.M with the given input arguments.
%
%      PHOEBE('Property','Value',...) creates a new PHOEBE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before phoebe_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to phoebe_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help phoebe

% Last Modified by GUIDE v2.5 23-Sep-2018 11:30:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @phoebe_OpeningFcn, ...
                   'gui_OutputFcn',  @phoebe_OutputFcn, ...
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


% --- Executes just before phoebe is made visible.
function phoebe_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to phoebe (see VARARGIN)

containing_dir = mfilename('fullpath');
cd(containing_dir(1:end-7))
% Add that folder plus all subfolders to the path.
addpath(pwd);
addpath(genpath([pwd filesep 'subfunctions']));

%% Choose default command line output for phoebe
handles.output = hObject;

%% Create Patriot and Phoebe parameters tabs
create_tabs;

%% Loads predefined head atlas
load_atlas;

%% Loads initialization parameters (last saved Phoebe settings)
load([pwd filesep 'subfunctions' filesep 'init_blank.mat']);
if exist([pwd filesep 'init.mat'],'file') == 2
    load([pwd filesep 'init.mat']);
    if isempty(dig_pts_path) || ~isstr(dig_pts_path) || ~exist(dig_pts_path, 'file')
      warningMessage = sprintf('Warning: file does not exist:\n%s', dig_pts_path);
      uiwait(msgbox(warningMessage));
      load([pwd filesep 'subfunctions' filesep 'init_blank.mat']);
    end
    if isempty(pairings_path) || ~isstr(pairings_path) || ~exist(pairings_path, 'file')
      warningMessage = sprintf('Warning: file does not exist:\n%s', pairings_path);
      uiwait(msgbox(warningMessage));
      load([pwd filesep 'subfunctions' filesep 'init_blank.mat']);
    end
end


%% Check Version
if w3conn()
    latest_ver = str2double(webread('http://polloninilab.com/version.txt'));
    if current_ver < latest_ver
        uiwait(msgbox(['PHOEBE Version ' num2double(latest_ver) ' is available. Please download it at http://bitbucket.org/lpollonini/phoebe'],'PHOEBE update','help'))
%        choice = questdlg('A new version of PHOEBE is available. Would you like to update it?','Settings','Yes','No','Yes');
%         if strcmp(choice,'Yes') % download new version of phoebe
%             current_ver = latest_ver;
%             save([pwd filesep 'init.mat'],'current_ver','-append');
%             !phoebe_update.exe
%             exit
%             % Interrupt this execution
%         end
    end
end

%% Check if this is the first use: if yes, ask for configuration to save in 'init.mat'
handles.opacity = opacity;
handles.zoom_index = zoom_index;
if first_time == 1
    % Patriot user?
    choice = questdlg('Are you using Polhemus Patriot to digitize optodes?','Settings','Yes','No','Yes');
    if strcmp(choice,'Yes') % Patriot user
        patriot_user = 1;
        save([pwd filesep 'init.mat'],'patriot_user','-append');  % Saving patriot flag in init.mat
        
        choice = questdlg({'Do you have a previous 3D digitization to be used as default?','NOTE: This is for first time use only. The default digitization can be changed at any time in menu Settings'},'Digitized Optodes','Yes','No','Yes');
        if strcmp(choice,'Yes')
            [FileName,PathName] = uigetfile('*.txt','Select the subject digitization file','C:\Users\owner\Desktop\Digitization\*.txt');
            dig_pts_path = [PathName FileName];
            save([pwd filesep 'init.mat'],'dig_pts_path','-append')
            handles = load_dig_pts(handles,dig_pts_path); % Import digitized layout and transform into atlas space
            
            choice = questdlg({'Do you have a source-detector pairings file (Homer format) to be used as default?','If not, it will be created for you based on the S-D euclidean distance.'},'Digitized Optodes','Yes','No','Yes');
            if strcmp(choice,'Yes')
                [FileName,PathName] = uigetfile('*.SD','Select the optodes pairings file','*.SD');
                load([PathName FileName],'-mat');
                % maybe a sanity check to see if digitization and SD pairs go along?
                pairings_path = [PathName FileName];
                save([pwd filesep 'init.mat'],'pairings_path','-append')
            else % Not a default SD pairing
                % Choose range
                % Execute SDrange
                % Save SD file in chosen location
            end
            % Loads SD pairs to be considered for Phoebe, and saves info into GUI handler
            handles = load_SD(handles,pairings_path);
            % In the settings dialog, we can ask to recompute based on new distance
            plot_atlas
        
        else
            uiwait(warndlg('Since there is no default optode layout, please digitize your headgear before monitoring the scalp coupling with Phoebe','No default digitization file','modal'));
            plot_atlas_empty
        end
        %maybe settings window based on init file, so user can choose to
        %default one 3D scan, point to a different one, or perform a new probe clicking.
    
    else    % Non-Patriot user
        patriot_user = 0;
        save([pwd filesep 'init.mat'],'patriot_user','-append');  % Saving patriot flag in init.mat
        
        choice = questdlg('Do you want to use a previously saved 3D digitization?','Digitized Optodes','Yes','No','Yes');
        if strcmp(choice,'Yes')
            [FileName,PathName] = uigetfile('*.txt','Select the subject digitization file','C:\Users\owner\Desktop\Digitization\*.txt');
            dig_pts_path = [PathName FileName];
            save([pwd filesep 'init.mat'],'dig_pts_path','-append')
            handles = load_dig_pts(handles,dig_pts_path); % Import digitized layout and transform into atlas space
            
            choice = questdlg({'Do you have a source-detector pairings file (Homer format)?','If not, it will be created for you based on the S-D euclidean distance'},'Digitized Optodes','Yes','No','Yes');
            if strcmp(choice,'Yes')
                [FileName,PathName] = uigetfile('*.SD','Select the optodes pairings file','*.SD');
                load([PathName FileName],'-mat');
                pairings_path = [PathName FileName];
                save([pwd filesep 'init.mat'],'pairings_path','-append')
            else % Not a default SD pairing
                det_pts = handles.det_pts;
                src_pts = handles.src_pts;
                idx_min_cell=rangesearch(det_pts,src_pts,min_pts_range);
                idx_max_cell=rangesearch(det_pts,src_pts,max_pts_range);
                for i = 1:size(idx_max_cell,1)
                    idx{i,1}=setdiff(idx_max_cell{i,1},idx_min_cell{i,1});
                end
                row=1;
                MeasList=[];
                for i=1:size(src_pts,1)
                    det_conn=sort(idx{i,1});
                    for j=1:length(det_conn)
                        MeasList(row,:)=[i det_conn(j) 1 1];
                        row=row+1;
                    end
                end
                % convert for NIRS toolbox format
                SD.SrcPos = src_pts;
                SD.DetPos = det_pts;
                SD.nSrcs = size(src_pts,1);
                SD.nDets = size(det_pts,1);
                MeasList2 = MeasList;
                MeasList2(:,4) = 2;
                SD.MeasList = [MeasList; MeasList2];
                pairings_path = [PathName FileName(1:end-4) '_pairings.SD'];
                save(pairings_path,'SD');
                save([pwd filesep 'init.mat'],'pairings_path','-append')
            end
            % Loads SD pairs to be considered for Phoebe, and saves info into GUI handler
            handles = load_SD(handles,pairings_path); % If there is a digitization file, there must be a pairings path, so no need ot check for it as in my Phoebe_dist?
            % In the settings dialog, we can ask to recompute based on new distance
            plot_atlas
        
        else % No default digitization: point and click on head model?           
            choice = questdlg('Do you want to create an aproximate 3D digitization by locating optodes on a head model with mouse clicks?','Digitized Optodes','Yes','No','Yes');
            if strcmp(choice,'Yes')
                uiwait(probe_configuration(handles))
            else
                uiwait(warndlg('Since there is no default optode layout, please select a default Patriot digitization or create an aproximate one before monitoring the scalp coupling with Phoebe','No default digitization file','modal'));
            end
            plot_atlas_empty
        end
    end
    first_time = 0;
    save([pwd filesep 'init.mat'],'first_time','dig_pts_path','pairings_path','-append');  % Update first time flag in init.mat
    
else % After first time
    
    % Prepare one or two head models for plotting
    if double_view==0
        set(handles.uipanel_head,'SelectedObject',handles.radiobutton_singleview);
    else
        set(handles.uipanel_head,'SelectedObject',handles.radiobutton_doubleview);
    end

    % Import digitized layout and transform into atlas space
    if ~strcmp(dig_pts_path,'')
        handles = load_dig_pts(handles,dig_pts_path);
        % Loads SD pairs to be considered for Phoebe, and saves info into GUI handler
        handles = load_SD(handles,pairings_path); % If there is a digitization file, there must be a pairings path, so no need ot check for it as in my Phoebe_dist?
        % In the settings dialog, we can ask to recompute based on new distance
        handles.dig_pts_path = dig_pts_path;
        handles.pairings_path = pairings_path;
        plot_atlas
    else
        plot_atlas_empty
    end
    
%     switch device
%         case 'NIRx'
%            handles.nirx_cfg_path = nirx_cfg_path; 
%     end
    
end % End first or non-first time 

if ~isfield(handles,'src_num')
    set(handles.source_dropdown,'Value',1);
    set(handles.detector_dropdown,'Value',1);
    set(handles.fiducial_dropdown,'Value',1);
else
    set(handles.source_dropdown,'Value',handles.src_num);
    set(handles.detector_dropdown,'Value',handles.det_num);
    if size(handles.fid_pts,1)==3
        set(handles.fiducial_dropdown,'Value',1);
    else
        set(handles.fiducial_dropdown,'Value',2);
    end
end

set(handles.edit_lowcutoff,'string',num2str(fcut_min));
set(handles.edit_highcutoff,'string',num2str(fcut_max));
set(handles.edit_threshold,'string',num2str(sci_threshold));
set(handles.min_optode_dist_edit,'String',num2str(min_pts_range));
set(handles.max_optode_dist_edit,'String',num2str(max_pts_range));
set(handles.edit_SCIwindow,'string',num2str(sci_window));
set(handles.edit_spectral_threshold,'string',num2str(power_threshold));
set(handles.slider_opacity,'value',opacity);
handles.baud_rate = baud_rate;
handles.com_port = com_port;
handles.patriot_hemisphere = patriot_hemisphere;
handles.current_ver = current_ver;

% Update GUI handler and launch GUI inferface
guidata(hObject, handles);
% UIWAIT makes phoebe wait for user response (see UIRESUME)
% uiwait(handles.figure_main);


% Updates graphics when toggling between single and double view
function uipanel_head_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_head 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if eventdata.NewValue == handles.radiobutton_singleview 
    cla(handles.axes_right)
end
if eventdata.NewValue == handles.radiobutton_doubleview
    cla(handles.axes_left)
    if isfield(handles,'fid_pts')
        plot_atlas
    else
        plot_atlas_empty
    end
end
hold(handles.axes_right,'off')
set(handles.togglebutton_scan,'String','START MONITOR');
set(handles.radiobutton_singleview,'Enable','on');
set(handles.radiobutton_doubleview,'Enable','on');
guidata(hObject, handles);


%% MENU PROBE -> CREATE
function menu_probe_create_Callback(hObject, ~, handles)
uiwait(probe_configuration(handles))
cla(handles.axes_left)
cla(handles.axes_right)
handles = guidata(hObject);
handles = load_dig_pts(handles,handles.dig_pts_path);  % Import digitized layout and transform into atlas space
handles = load_SD(handles,handles.pairings_path);  % Loads SD pairs to be considered for Phoebe, and saves info into GUI handler
plot_atlas


%% START/STOP MONITORING
function togglebutton_scan_Callback(hObject, ~, handles, FileName, PathName)
% Hint: get(hObject,'Value') returns toggle state of togglebutton_scan

if get(hObject,'Value') %If currently STOP
    set(handles.togglebutton_scan,'String','STOP MONITOR');
    set(handles.radiobutton_singleview,'Enable','off');
    set(handles.radiobutton_doubleview,'Enable','off');
    drawnow
    
    lib = lsl_loadlib(); % Loading LSL library
    if isempty(lib)
       uiwait(msgbox('LSL library not loaded','PHOEBE','error')) 
    end
    
    result = lsl_resolve_byprop(lib,'name','NIRStar');
    if ~isempty(result)
        inlet = lsl_inlet(result{1});
        [dummy,~] = inlet.pull_chunk();
    else
        uiwait(warndlg('Please PREVIEW or RECORD data in NIRStar and ensure that the LSL streaming (full mask) is active','PHOEBE'))
        set(handles.togglebutton_scan,'String','START MONITOR');
        set(handles.togglebutton_scan,'Value',0);
        set(handles.radiobutton_singleview,'Enable','on');
        set(handles.radiobutton_doubleview,'Enable','on');
        guidata(hObject,handles)
        return
    end
    
    switch(get(handles.popupmenu_device,'value'))   % Instrument
        case 1  % NIRx
            if (result{1}.channel_count - 1) ~= 2*size(handles.src_pts,1)*size(handles.det_pts,1)
                uiwait(warndlg('The size of streaming data does not correspond to the probe size. Please ensure that the number of sources and detectors in NIRStar and PHOEBE are set consistently.','PHOEBE'))
                set(handles.togglebutton_scan,'String','START MONITOR');
                set(handles.togglebutton_scan,'Value',0);
                set(handles.radiobutton_singleview,'Enable','on');
                set(handles.radiobutton_doubleview,'Enable','on');
                guidata(hObject,handles)
                return
            end
            handles.fs = result{1}.nominal_srate;
            buffer_rows = floor(handles.fs*str2double(get(handles.edit_SCIwindow,'string'))); % Number of frames to be buffered
            lsl_buffer = zeros(buffer_rows,result{1}.channel_count);
        case 2  % Another device
            
    end
    num_optodes = size(handles.src_pts,1) + size(handles.det_pts,1); % Computes total optodes based on loaded probe 
    
    % Plots filled markers, but we need the handlers h1 and h2 of optodes to update their colors below (weird)
    hold(handles.axes_left,'on')
    hold(handles.axes_right,'on')
    delete(handles.h_src_left);
    delete(handles.h_det_left);
    h1=scatter3(handles.axes_left,handles.src_pts(:,1),handles.src_pts(:,2),handles.src_pts(:,3),60,'r','fill');
    h2=scatter3(handles.axes_left,handles.det_pts(:,1),handles.det_pts(:,2),handles.det_pts(:,3),60,'b','fill','s');
    if get(handles.uipanel_head,'SelectedObject') == handles.radiobutton_doubleview
        delete(handles.h_src_right);
        delete(handles.h_det_right);
        h3=scatter3(handles.axes_right,handles.src_pts(:,1),handles.src_pts(:,2),handles.src_pts(:,3),60,'r','fill');
        h4=scatter3(handles.axes_right,handles.det_pts(:,1),handles.det_pts(:,2),handles.det_pts(:,3),60,'b','fill','s'); 
    end
    
    % Reads the filter parameters form panel and computes low-pass coeffs
    % (FIX make script or function)
    fcut_min=str2double(get(handles.edit_lowcutoff,'string'));
    fcut_max=str2double(get(handles.edit_highcutoff,'string'));
    if fcut_max >= (handles.fs)/2
        fcut_max = (handles.fs)/2 - eps;
        set(handles.edit_highcutoff,'string',num2str(fcut_max));
        uiwait(warndlg('The highpass cutoff has been reduced to Nyquist sampling rate. This setting will be saved for future use.'));
    end
    [B1,A1]=butter(1,[fcut_min*(2/handles.fs) fcut_max*(2/handles.fs)]);

else % Stops monitoring
    set(handles.togglebutton_scan,'String','START MONITOR');
    set(handles.radiobutton_singleview,'Enable','on');
    set(handles.radiobutton_doubleview,'Enable','on');
end

% If START/STOP button is START mode, run this loop indefinitely until the button is toggled
while ishandle(hObject) && get(hObject,'Value')
    switch(get(handles.popupmenu_device,'value'))
        case 1  % NIRx
            while nnz(lsl_buffer(:,1)) < size(lsl_buffer,1) % if the FIFO buffer is not totally full (first few secs), fill it
                while 1
                    [chunk,~] = inlet.pull_chunk(); % Pull a chunk of fresh samples
                    if ~isempty(chunk)
                        break
                    end
                end
                chunk = chunk';
                if size(chunk,1) < size(lsl_buffer,1)    % If the data chunk is smaller than the buffer
                    lsl_buffer(1:size(lsl_buffer,1)-size(chunk,1),:) = lsl_buffer(size(chunk,1)+1:end,:); % Shift up the buffer to make room 
                    lsl_buffer(size(lsl_buffer,1)-size(chunk,1)+1:end,:) = chunk; % Put chunk in buffer
                else
                    lsl_buffer(:,:) = chunk(size(chunk,1)-size(lsl_buffer,1)+1:end,:); % Put chunk in buffer
                end% Put chunk in buffer
            end
            
            while 1
                [chunk,~] = inlet.pull_chunk(); % Pull a chunk of fresh samples
                if ~isempty(chunk)
                    break
                end
            end
            chunk = chunk';
            if size(chunk,1) < size(lsl_buffer,1)    % If the data chunk is smaller than the buffer
                lsl_buffer(1:size(lsl_buffer,1)-size(chunk,1),:) = lsl_buffer(size(chunk,1)+1:end,:); % Shift up the buffer to make room 
                lsl_buffer(size(lsl_buffer,1)-size(chunk,1)+1:end,:) = chunk; % Put chunk in buffer
            else
                lsl_buffer(:,:) = chunk(size(chunk,1)-size(lsl_buffer,1)+1:end,:); % Put chunk in buffer
            end
            
            nirs_data1 = lsl_buffer(:,2:size(handles.src_pts,1)*size(handles.det_pts,1)+1);
            nirs_data2 = lsl_buffer(:,size(handles.src_pts,1)*size(handles.det_pts,1)+2:end);
            
        case 2  % Other device
    
    end % End readout of incoming data
    
    % Filter everything bu the cardiac component
    filtered_nirs_data1=filtfilt(B1,A1,nirs_data1);       % Cardiac bandwidth
    filtered_nirs_data1=filtered_nirs_data1./repmat(std(filtered_nirs_data1,0,1),size(filtered_nirs_data1,1),1); % Normalized heartbeat
    filtered_nirs_data2=filtfilt(B1,A1,nirs_data2);       % Cardiac bandwidth
    filtered_nirs_data2=filtered_nirs_data2./repmat(std(filtered_nirs_data2,0,1),size(filtered_nirs_data2,1),1); % Normalized heartbeat
    
    % Distributes all the measures in the optode "battlefield" matrix
    sci_matrix = zeros(num_optodes,num_optodes);    % Number of optode is from the user's layout, not the machine
    power_matrix = zeros(num_optodes,num_optodes);
    fpower_matrix = zeros(num_optodes,num_optodes);
    A = zeros(num_optodes,num_optodes);
    for i=1:size(handles.SDpairs,1)
        col = (handles.SDpairs(i,1)-1)*handles.det_num + handles.SDpairs(i,2);  
        similarity = xcorr(filtered_nirs_data1(:,col),filtered_nirs_data2(:,col),'unbiased');  %cross-correlate the two wavelength signals - both should have cardiac pulsations
        similarity = length(filtered_nirs_data1(:,col))*similarity./sqrt(sum(abs(filtered_nirs_data1(:,col)).^2)*sum(abs(filtered_nirs_data2(:,col)).^2));  % this makes the SCI=1 at lag zero when x1=x2 AND makes the power estimate independent of signal length, amplitude and Fs
        [pxx,f] = periodogram(similarity,hamming(length(similarity)),length(similarity),handles.fs,'power');
        [pwrest,idx] = max(pxx(f<1.7)); % FIX Make it age-dependent
        sci=similarity(length(filtered_nirs_data1(1:end,i)));
        power=pwrest;
        fpower=f(idx);
        sci_matrix(handles.SDpairs(i,1),handles.src_num+handles.SDpairs(i,2)) = sci;    % Adjust not based on machine
        power_matrix(handles.SDpairs(i,1),handles.src_num+handles.SDpairs(i,2)) = power;
        fpower_matrix(handles.SDpairs(i,1),handles.src_num+handles.SDpairs(i,2)) = fpower;
        A(handles.SDpairs(i,1),handles.src_num+handles.SDpairs(i,2)) = 1;   % Adjacency matrix: marks all the active pairs with 1, rest is 0
    end
    W = (sci_matrix >= str2double(get(handles.edit_threshold,'string'))) & (power_matrix >= str2double(get(handles.edit_spectral_threshold,'string'))); % Weight boolean matrix: here we set the criteria for pass
    % Computes optodes coupling status: coupled (1), uncoupled (0) or undetermined (-1).
    [optodes_status] = boolean_system(num_optodes,A,W); 
    optodes_color = zeros(length(optodes_status),3);
    for i=1:length(optodes_status)
        switch(optodes_status(i))
            case 1
                optodes_color(i,:) = [0 1 0];
            case 0
                optodes_color(i,:) = [1 0 0];
            case -1
                optodes_color(i,:) = [1 1 0];
        end
    end
   
    set(h1,'CData',optodes_color(1:handles.src_num,:));
    set(h2,'CData',optodes_color(handles.src_num+1:end,:));
    if get(handles.uipanel_head,'SelectedObject')==handles.radiobutton_doubleview
        set(h3,'CData',optodes_color(1:handles.src_num,:));
        set(h4,'CData',optodes_color(handles.src_num+1:end,:));
    end
    drawnow
end


%% CLOSE FUNCTION
% When user attempts to close GUI, then save all the parameters from the panel
function figure_main_CloseRequestFcn(hObject, ~, handles)

fcut_min = str2double(get(handles.edit_lowcutoff,'string'));
fcut_max = str2double(get(handles.edit_highcutoff,'string'));
sci_threshold = str2double(get(handles.edit_threshold,'string'));
if get(handles.uipanel_head,'SelectedObject')==handles.radiobutton_singleview
    double_view = 0;
else
    double_view = 1;
end
sci_window = str2double(get(handles.edit_SCIwindow,'string'));
power_threshold = str2double(get(handles.edit_spectral_threshold,'string'));
%com_port = get(handles.com_port_edit,'String');
%baud_rate = get(handles.baud_rate_edit,'String');
min_pts_range = str2double(get(handles.min_optode_dist_edit,'String'));
max_pts_range = str2double(get(handles.max_optode_dist_edit,'String'));
opacity = get(handles.slider_opacity,'Value');
save([pwd filesep 'init.mat'],'fcut_min','fcut_max','sci_threshold','double_view','sci_window','sci_window','power_threshold','min_pts_range','max_pts_range','opacity','-append')
% Hint: delete(hObject) closes the figure
if ~isempty(instrfind)
    fclose(instrfind)
end
set(handles.togglebutton_scan,'String','START MONITOR');
set(handles.radiobutton_singleview,'Enable','on');
set(handles.radiobutton_doubleview,'Enable','on');
delete(hObject);


%% SLIDER OPACITY
function slider_opacity_Callback(hObject, ~, handles)
% hObject    handle to slider_opacity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.opacity = get(hObject,'Value');
cla(handles.axes_left);
cla(handles.axes_right);
if isfield(handles,'fid_pts')
    plot_atlas;
else
    plot_atlas_empty;
end
guidata(hObject, handles);


%% SLIDER ZOOM
% --- Executes on slider movement.
function slider_zoom_Callback(hObject, eventdata, handles)
% hObject    handle to slider_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.zoom_index = 1 - get(hObject,'Value');
cla(handles.axes_left);
cla(handles.axes_right);
if isfield(handles,'fid_pts')
    plot_atlas;
else
    plot_atlas_empty;
end
guidata(hObject, handles);


%% LOAD PROBE
function menu_probe_load_Callback(hObject, ~, handles)
%
set(handles.uipanel_head,'SelectedObject',handles.radiobutton_doubleview);
cla(handles.axes_left);
cla(handles.axes_right);
[FileName,PathName] = uigetfile('*.txt','Select the subject digitization file','C:\Users\owner\Desktop\Digitization\*.txt');
dig_pts_path = [PathName FileName];
save([pwd filesep 'init.mat'],'dig_pts_path','-append')
handles.dig_pts_path = dig_pts_path;
choice = questdlg({'Do you have a source-detector pairings file (Homer format) to be used as default?','If not, it will be created for you based on the S-D euclidean distance.'},'Digitized Optodes','Yes','No','Yes');
if strcmp(choice,'Yes')
    [FileName,PathName] = uigetfile('*.SD','Select the optodes pairings file',[PathName '*.SD']);
    load([PathName FileName],'-mat');
    % maybe a sanity check to see if digitization and SD pairs go along?
    pairings_path = [PathName FileName];
    save([pwd filesep 'init.mat'],'pairings_path','-append')
    handles.pairings_path = pairings_path;
else
    guidata(hObject, handles);
    uiwait(create_pairings(handles))
    handles_main = findobj('Tag','figure_main');
    handles = guidata(handles_main);
end
handles = load_dig_pts(handles,handles.dig_pts_path);  % Import digitized layout and transform into atlas space
handles = load_SD(handles,handles.pairings_path);  % Loads SD pairs to be considered for Phoebe, and saves info into GUI handler
plot_atlas
guidata(hObject, handles);


%% MENU TOOLS -> CONVERT PATRIOT TO MNI
% --------------------------------------------------------------------
function menu_patriot2mni_Callback(hObject, ~, handles)
% hObject    handle to menu_patriot2mni (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'dig_pts_path')
    current_dig_pts_path = handles.dig_pts_path;
    [current_dig_pts_path,~,~] = fileparts(current_dig_pts_path);
else
    current_dig_pts_path = [pwd filesep '*.txt'];
end
[FileName,PathName] = uigetfile('*.txt','Select the digitization file',current_dig_pts_path);
dig_pts_path = [PathName FileName];
handles = load_dig_pts(handles,dig_pts_path);  % Import digitized layout and transform into atlas space

[FileName,PathName] = uiputfile('*.txt','Saving the subject digitization file',[PathName FileName(1:end-4) '_mni.txt']);
fileID = fopen([PathName FileName],'w');
fprintf(fileID,'nz: %.1f %.1f %.1f\n',handles.fid_pts(1,1),handles.fid_pts(1,2),handles.fid_pts(1,3));
fprintf(fileID,'ar: %.1f %.1f %.1f\n',handles.fid_pts(2,1),handles.fid_pts(2,2),handles.fid_pts(2,3));
fprintf(fileID,'al: %.1f %.1f %.1f\n',handles.fid_pts(3,1),handles.fid_pts(3,2),handles.fid_pts(3,3));

for i = 1:length(handles.src_pts)
    fprintf(fileID,'s%d: %.1f %.1f %.1f\n',i,handles.src_pts(i,1),handles.src_pts(i,2),handles.src_pts(i,3));
end
for i = 1:1:length(handles.det_pts)
    fprintf(fileID,'d%d: %.1f %.1f %.1f\n',i,handles.det_pts(i,1),handles.det_pts(i,2),handles.det_pts(i,3));
end
fclose(fileID);


%% MENU TOOLS -> NEW RANGE
% --------------------------------------------------------------------
function menu_range2pairings_Callback(hObject, ~, handles)
% hObject    handle to menu_range2pairings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





%% MENU ABOUT -> CHECK UPDATE
function menu_check_update_Callback(hObject, ~, handles)
% hObject    handle to menu_check_update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
current_ver = handles.current_ver;
h_msgchk = msgbox('Checking for updated version...','PHOEBE');
latest_ver = str2double(webread('https://www.uh.edu/tech/pollonini/phoebe/version.txt'));
if current_ver < latest_ver
    close(h_msgchk);
    uiwait(msgbox('A new version of PHOEBE is available on Bitbucket','PHOEBE','replace'));
%     choice = questdlg('A new version of PHOEBE is available. Would you like to update it?','Settings','Yes','No','Yes');
%     if strcmp(choice,'Yes') % download new version of phoebe
%         current_ver = latest_ver;
%         save([pwd filesep 'init.mat'],'current_ver','-append');
%         !phoebe_update.exe&
%         exit
%         %Either exit the application from here, or let Update kill it
%     end
else
    uiwait(msgbox('PHOEBE is up to date','PHOEBE','replace'));
end

%% DIGITIZE BUTTON
% --- Executes on button press in pushbutton_digitize.
function [handles] = pushbutton_digitize_Callback(hObject, ~, handles)
% hObject    handle to pushbutton_digitize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% patriotGUIcodeREALTIME

[handles] = patriot(handles);
guidata(hObject,handles); 


%% MENU PREFERENCES
% --------------------------------------------------------------------
function menu_settings_Callback(hObject, ~, handles)
% hObject    handle to menu_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(phoebe_settings(handles))



%% GARBAGE GUI FUNCTIONS

% Outputs from this function are returned to the command line.
function varargout = phoebe_OutputFcn(hObject, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function menu_online_help_Callback(hObject, ~, handles)
% hObject    handle to menu_online_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
url = 'https://bitbucket.org/lpollonini/phoebe/wiki';
web(url,'-browser')

% --------------------------------------------------------------------
function menu_about_Callback(hObject, ~, handles)
% hObject    handle to menu_about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_tools_Callback(hObject, eventdata, handles)
% hObject    handle to menu_tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_probe_Callback(hObject, eventdata, handles)
% hObject    handle to menu_probe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_help_Callback(hObject, eventdata, handles)
% hObject    handle to menu_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_digitizer_Callback(hObject, eventdata, handles)
% hObject    handle to menu_digitizer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function varargout = SpectralGUI_BackUp(varargin)
% SPECTRALGUI MATLAB code for SpectralGUI.fig
%      SPECTRALGUI, by itself, creates a new SPECTRALGUI or raises the existing
%      singleton*.
%
%      H = SPECTRALGUI returns the handle to a new SPECTRALGUI or the handle to
%      the existing singleton*.
%
%      SPECTRALGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPECTRALGUI.M with the given input arguments.
%
%      SPECTRALGUI('Property','Value',...) creates a new SPECTRALGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SpectralGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SpectralGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SpectralGUI

% Last Modified by GUIDE v2.5 08-Sep-2017 08:43:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SpectralGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SpectralGUI_OutputFcn, ...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes just before SpectralGUI is made visible.
function SpectralGUI_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SpectralGUI (see VARARGIN)

% Choose default command line output for SpectralGUI
handles.output = hObject;

%%%%%%%%%%%
%Initialze variables here 
cla(handles.Graph_1); cla(handles.Graph_2); cla(handles.Graph_3); 
cla(handles.axes5); cla(handles.axes6);

handles = initializing(hObject, handles);

set(handles.get_status, 'String','Spectral Clustering by Tram Pham, UC Berkeley ');

handles.DataExtensions = {'*.m;*.mlx;*.fig;*.mat;*.slx;*.mdl; *.csv',...
 'MATLAB Files (*.m,*.mlx,*.fig,*.mat,*.slx,*.mdl)';
   '*.m;*.mlx',  'Code files (*.m,*.mlx)'; ...
   '*.fig','Figures (*.fig)'; ...
   '*.mat','MAT-files (*.mat)'; ...
   '*.mdl;*.slx','Models (*.slx, *.mdl)'; ...
   '*.*',  'All Files (*.*)';
   '*.csv', 'CSV FILEs'};

set(handles.Normalize, 'Enable', 'off')

handles.fileName = '';
handles.command=[];
handles.pathname = '';
handles.currentDataSet = 0;
handles.current_eps = 1;
handles.current_sig = 1;
handles.numOfNeighbors = 2;
handles.clusterType = 2;
handles.L = []; handles.U = []; handles.C = [];
handles.eigenvalues = [];
handles.W = [];
handles.chosen_dimensions = 1:2;
handles.DATA = [];

global Spectral;


Spectral = struct('CurrentNumPoints',150, 'MinNumPoints',10, 'MaxNumPoints',500, ...
                   'CurrentDataDim',2, 'MinDataDim',2, 'MaxDataDim',200,... 
                   'CurrentSimFct','Gaussian kernel',...
                   'CurrentGraph1', 1, ... %refers to position in list, 1=eps
                   'CurrentGraphPara1',0.5,... 
                   'DefaultK',5,'DefaultEps',0.5,... %default para values
                   'CurrentSigma',0.5, 'MinLogSigma',-2, 'MaxLogSigma',2, ... %sigma of the similarity function
                   'DefaultSigma',0.5, ... 
                   'CurrentVarNoise',0.01, 'MinVarNoise',0,'MaxVarNoise',0.1,... 
                    'MinK', 1, 'MaxK',50, 'MinEps', 0, 'MaxEps',100, ...
                   'points', [], 'labels', [] , ... %data set
                   'Distane', [], ... % data distance matrix
                   'Similarity', [], ...
                   'W', []); %data similarity matrix
                   
               
               
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SpectralGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
% function varargout = SpectralGUI_OutputFcn(hObject, eventdata, handles) 
function varargout = SpectralGUI_OutputFcn(~, ~, handles)
varargout{1} = handles.output;

% --- Executes on button press in LoadData.
function LoadData_Callback(hObject, ~, handles)
% hObject    handle to LoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.select_sampleData, 'Enable', 'off');


[handles.fileName, handles.pathname] =  uigetfile(handles.DataExtensions,...
    'Select a MATLAB .mat file', '');

if ~isequal(handles.fileName,0)
   set(handles.get_status, 'String', 'Busy Loading Data');
   disp(['User selected ', fullfile(handles.pathname, handles.fileName)])
   
% handles.command = sprintf('load(''%s'')', handles.fileName);
% handles.DATA = load(handles.fileName);
% handles.field = fieldnames(handles.DATA);
% name = handles.field{1};
% handles.DATA = handles.DATA(name);
handles.DATA = csvread(fullfile(relativepath(handles.pathname), ...
                        handles.fileName));
 
   set(handles.get_status, 'String', 'Done Loading Data')
   
   [~, handles.dataName, handles.DataExtension] = fileparts(handles.fileName);
   
   %update name of data and its info.
   set(handles.DataName, 'String', handles.dataName);
   set(handles.textDataSize, 'String', int2str(size(handles.DATA,1)));
   set(handles.textDimension, 'String', int2str(size(handles.DATA,2)));
   
   set(handles.Normalize, 'Enable', 'on');
   set(handles.LoadData, 'Enable', 'off');
   
end
guidata(hObject, handles);



function Normalize_Callback(hObject, eventdata, handles)
% hObject    handle to Normalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.get_status, 'String', 'Busy Normalizing Data Set');
handles.DATA = normalizeDATA(handles.DATA);

set(handles.get_status, 'String', 'Done Normalizaing Data Set');

guidata(hObject, handles);

% --- Executes on selection change in select_sampleData.
function select_sampleData_Callback(hObject, eventdata, handles)
set(handles.LoadData, 'Enable', 'off');
global Spectral;

% global Spectral;
handles.currentDataSet = get(hObject, 'Value');
if handles.currentDataSet == 1
   set(handles.DataName, 'String', 'Two Moons');
   set(handles.textDataSize, 'String', int2str(Spectral.CurrentNumPoints));
   set(handles.textDimension, 'String', int2str(Spectral.CurrentDataDim));
elseif handles.currentDataSet == 2
   set(handles.DataName, 'String', 'Two Gaussians');
   set(handles.textDataSize, 'String', int2str(Spectral.CurrentNumPoints));
   set(handles.textDimension, 'String', int2str(Spectral.CurrentDataDim));
end 
guidata(hObject, handles);


% --- Executes on selection change in popupmenu1_graphtype.
function popupmenu1_graphtype_Callback(hObject, eventdata, handles)
handles.currentGraph = get(hObject,'Value');

switch handles.currentGraph
    case 1 
        set(handles.numNeighbors, 'Enable', 'off');
    case 2 
        set(handles.current_epsilon, 'Enable', 'off');
        set(handles.numNeighbors, 'Enable', 'on');
    case 3 
        set(handles.current_epsilon, 'Enable', 'off');
        set(handles.numNeighbors, 'Enable', 'on');
    case 4 
        set(handles.current_epsilon, 'Enable', 'on');
        set(handles.numNeighbors, 'Enable', 'off');

    case 5
        set(handles.current_sigma, 'Enable', 'off');
        set(handles.numNeighbors, 'Enable', 'off');
        set(handles.current_epsilon, 'Enable', 'on');
    otherwise
end 

guidata(hObject, handles)



function numNeighbors_Callback(hObject, eventdata, handles)
handles.numOfNeighbors = str2double(get(hObject,'String'));

guidata(hObject,handles);



function current_epsilon_Callback(hObject, eventdata, handles)
% hObject    handle to current_epsilon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of current_epsilon as text
%        str2double(get(hObject,'String')) returns contents of current_epsilon as a double
handles.current_eps = str2double(get(hObject,'String'));

guidata(hObject,handles);


function current_sigma_Callback(hObject, eventdata, handles)
global Spectral;
handles.current_sig = str2double(get(hObject,'String'));

Spectral.CurrentSigma = handles.current_sig;

guidata(hObject,handles);




% --- Executes on selection change in popupmenu_cluster_type.
function popupmenu_cluster_type_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_cluster_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_cluster_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_cluster_type
handles.clusterType = get(hObject, 'Value');

guidata(hObject, handles);


function NumCluster_Callback(hObject, eventdata, handles)
% hObject    handle to NumCluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumCluster as text
%        str2double(get(hObject,'String')) returns contents of NumCluster as a double
handles.current_NumCluster = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes on button press in graph_button.
function graph_button_Callback(hObject, eventdata, handles)
% hObject    handle to graph_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Plot_Button_Graph(handles)



% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



%%%%%%%%%%%%%%%%%%%%%%% Helpers Functions %%%%%%%%%%%%%
function Plot_Button_Graph(handles)

global Spectral;

%clear all graph 
cla(handles.Graph_1); cla(handles.Graph_2); cla(handles.Graph_3); 
cla(handles.axes5); cla(handles.axes6);
drawnow

if ~isequal(handles.currentDataSet,0)
     
switch handles.currentDataSet
  
  case 1 % Two Moons with balanced classes [0.5,0.5]
  density = 2;
  [pointstrans,~] =GD_GenerateData(density,Spectral.CurrentNumPoints, Spectral.CurrentDataDim,...
                                [0.5,0.5], Spectral.CurrentVarNoise);
                                 
 
 case 2 % Two isotropic Gaussians balanced classes [0.5,0.5]
  density=3;
  [pointstrans,~] =GD_GenerateData(density,Spectral.CurrentNumPoints,...
                              Spectral.CurrentDataDim , [0.5,0.5,0],...                                                
                                             Spectral.CurrentVarNoise);
  

   
end

Spectral.points = pointstrans';

set(handles.get_status, 'String', 'Plotting Data');
set(gcf, 'CurrentAxes', handles.Graph_1);
plot_Data2D(Spectral.points);
title('Data Set First 2 Dimension');
set(handles.get_status, 'String', 'Done Plotting Data');

% second plot for connectivity 
Spectral.Distance = squaredDistance(Spectral.points, Spectral.points);
Spectral.Similarity = exp(- Spectral.Distance/ (2*Spectral.CurrentSigma^2));

%%% check current graph type to calculate adjancency matrix and plot it 

switch handles.currentGraph 
    case 1 %full similarity 
        Spectral.W = full_SimGraph(Spectral.points, handles.current_sig);
    case 2 %normal
        Spectral.W = kNN_SimGraph(Spectral.Similarity, handles.numOfNeighbors,...
                                    'normal', 1, handles.current_sig);
    case 3 %mutual 
        Spectral.W = kNN_SimGraph(Spectral.Similarity, handles.numOfNeighbors,...
                                    'mutual', 1, handles.current_sig);
    case 4 %weighted epsilon
        Spectral.W = epsilon_SimGraph(Spectral.Distance, handles.current_eps, 1);
        Spectral.W = Spectral.W *Spectral.Similarity;
    case 5  %unweighted epsilon
       Spectral.W = epsilon_SimGraph(Spectral.Distance, handles.current_eps, 2);
end
   
set(handles.get_status, 'String', 'Plotting Adjancency Matrix');
set(gcf,'CurrentAxes', handles.Graph_2); cla
imagesc(Spectral.W)
title('Adjancency Matrix');
set(handles.get_status, 'String', 'Done Plotting Adjancency Matrix');

% thrid plot for similarity graph 
set(handles.get_status, 'String', 'Plotting Similarity Matrix');
set(gcf, 'CurrentAxes', handles.Graph_3); 
GD_PlotGraph(Spectral.points(:, 1:2), Spectral.W, 'Similarity Graph');
set(handles.get_status, 'String', 'Done Plotting Similarity Matrix');


else

set(handles.get_status, 'String', 'Plotting Data');

set(gcf, 'CurrentAxes', handles.Graph_1);
% figure(2)
plot_Data2D(handles.DATA)
% plotData3D(handles.DATA);
title('Data Set First 2 Dimension');
set(handles.get_status, 'String', 'Done Plotting Data');

handles.DATA_holder = handles.DATA;
handles.distance_matrix = distance(handles.DATA_holder, handles.DATA_holder);
handles.similarity_matrix = Gaussian(handles.distance_matrix,handles.current_sig);

switch handles.currentGraph 
    case 1 %full similarity 
       Spectral.W = full_SimGraph(handles.DATA, handles.current_sig);
    case 2 %normal 
        Spectral.W = kNN_SimGraph(handles.similarity_matrix, handles.numOfNeighbors,...
                                    'normal', 1, handles.current_sig);
    case 3 %mutual 
        Spectral.W = kNN_SimGraph(handles.similarity_matrix, handles.numOfNeighbors,...
                                    'mutual', 1, handles.current_sig);
    case 4 %weighted epsilon
        Spectral.W = epsilon_SimGraph(handles.distance_matrix, handles.current_eps, 1);
        Spectral.W = Spectral.W*handles.similarity_matrix;
    case 5 %unweighted epsilon
        Spectral.W = epsilon_SimGraph(handles.distance_matrix, handles.current_eps, 2);
end

set(handles.get_status, 'String', 'Plotting Adjancency Matrix');

set(gcf,'CurrentAxes', handles.Graph_2); cla
imagesc(Spectral.W)
title('Adjancency Matrix');
set(handles.get_status, 'String', 'Done Plotting Adjancency Matrix');

set(handles.get_status, 'String', 'Plotting Similarity Matrix');

set(gcf, 'CurrentAxes', handles.Graph_3); 
GD_PlotGraph(handles.DATA(:, 1:2), Spectral.W, 'Similarity Graph');

set(handles.get_status, 'String', 'Done Plotting Similarity Matrix');

end

set(handles.select_sampleData, 'Enable', 'on');
set(handles.LoadData, 'Enable', 'on');
set(handles.Normalize, 'Enable', 'on');


 
% --- Executes on button press in Cluster_Data_Button.
function Cluster_Data_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Cluster_Data_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Spectral;
if size(handles.DATA,1) == 0
    handles.DATA = Spectral.points;
end 

[handles.C, handles.L, handles.U, handles.eigenvalues, handles.cidx] = Spectral_test(Spectral.W,...
                             handles.current_NumCluster,handles.clusterType );
                         
% handles.C = convertClusterVector(handles.C);                    
set(handles.get_status, 'String', 'Done Clustering Data');
guidata(hObject, handles);


                                    

% --- Executes on button press in PlotCluster_button.
function PlotCluster_button_Callback(hObject, eventdata, handles)
% hObject    handle to PlotCluster_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes during object creation, after setting all properties.
set(gcf, 'CurrentAxes', handles.axes5);
% plotClusteredData(handles.DATA, full(handles.C), handles.NumCluster)
% figure(1)
plotClusterStarCoordinates(hObject, handles);
% plotClusteredData(handles.DATA, full(handles.C), handles.current_NumCluster);
title('Clustered Data');




function get_dimensions_Callback(hObject, eventdata, handles)
% hObject    handle to get_dimensions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of get_dimensions as text
%        str2double(get(hObject,'String')) returns contents of get_dimensions as a double
handles.chosen_dimensions = get(hObject,'Value');


% --- Executes on button press in plot_disimilarity.
function plot_disimilarity_Callback(hObject, eventdata, handles)
% hObject    handle to plot_disimilarity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(gcf, 'CurrentAxes', handles.axes6);
silhouette(handles.DATA,handles.cidx)
title('Silhouette value for Dismilarity');
% scatter(handles.eigenvalues, 'bx');

% --- Executes on button press in tsNE.
function tsNE_Callback(hObject, eventdata, handles)
% hObject    handle to tsNE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Save_cluster_result.
function Save_cluster_result_Callback(hObject, eventdata, handles)
% hObject    handle to Save_cluster_result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





%%%%%%%%%%%%%%%%%%%%%%%% Empty Functions %%%%%%%%%%%%%%%%%%%%%%

function NumCluster_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumCluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object creation, after setting all properties.
function select_sampleData_CreateFcn(hObject, eventdata, handles)
% hObject    handle to select_sampleData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Normalize.


% --- Executes during object creation, after setting all properties.
function Normalize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Normalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function current_sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to current_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function popupmenu1_graphtype_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function numNeighbors_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numNeighbors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function popupmenu_cluster_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_cluster_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function current_epsilon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to current_epsilon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Save_clusteredData.
function Save_clusteredData_Callback(hObject, eventdata, handles)
% hObject    handle to Save_clusteredData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in Save_Data_Plot.
function Save_Data_Plot_Callback(hObject, eventdata, handles)
% hObject    handle to Save_Data_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Save_Adjancency.
function Save_Adjancency_Callback(hObject, eventdata, handles)
% hObject    handle to Save_Adjancency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Save_SimGraph.
function Save_SimGraph_Callback(hObject, eventdata, handles)
% hObject    handle to Save_SimGraph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function get_dimensions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to get_dimensions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function varargout = postprocessing(varargin)
% POSTPROCESSING M-file for postprocessing.fig
%      POSTPROCESSING, by itself, creates a new POSTPROCESSING or raises the existing
%      singleton*.
%
%      H = POSTPROCESSING returns the handle to a new POSTPROCESSING or the handle to
%      the existing singleton*.
%
%      POSTPROCESSING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POSTPROCESSING.M with the given input arguments.
%
%      POSTPROCESSING('Property','Value',...) creates a new POSTPROCESSING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before postprocessing_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to postprocessing_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help postprocessing

% Last Modified by GUIDE v2.5 13-Nov-2017 15:22:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @postprocessing_OpeningFcn, ...
                   'gui_OutputFcn',  @postprocessing_OutputFcn, ...
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

% --- Executes just before postprocessing is made visible.
function postprocessing_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to postprocessing (see VARARGIN)

% Choose default command line output for postprocessing
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes postprocessing wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = postprocessing_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% LOAD RESULTS
% --- Executes on button press in resButton.
function resButton_Callback(hObject, eventdata, handles)
% hObject    handle to resButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Choose the case file to run
[filename, pathname, ~] = uigetfile('*.mat', 'Choose results M-file');

% Load first structure in mat file as results
a = load([pathname , filename]);
b = fields(a);
results = a.(b{1});
handles.results = results;

% Now plot results in pareto front
plot(handles.axesPareto, results.fval(:,1) , results.fval(:,2) , '.b')

% Now sort results according to second objective
[handles.fval_sorted , handles.x_sorted ] = sort_pareto_data(results.fval , results.x);
% Flip to start from the thin airfoils
handles.fval_sorted = flipud(handles.fval_sorted);
handles.x_sorted    = flipud(handles.x_sorted);

% Now convert to Full space if optimization was ran in a subspace
if handles.GM.use_eqconstraint_reducer == true
    handles.x_sorted = handles.GM.ECR.array_from_null_to_full_space(handles.x_sorted);
end
% Otherwise do nothing!

% Make Cell matrix to fill Uitable
% a = mat2cell(handles.fval_sorted, ones(size(handles.fval_sorted,1),1), ones(1,size(handles.fval_sorted,2)));

% Make Cell array to fill Uitable (careful in keeping logicals as such, so
% do not convert them to float before going to cell)
uicell  = [num2cell(false(size(handles.fval_sorted, 1), 2)) , num2cell(handles.fval_sorted)];
% Write data into array
set(handles.uitable1, 'Data' , uicell);

% Set Cl bounds to start with
handles.mincl = -0.1;
handles.maxcl = 1.8;

handles.minalpha = 0;
handles.maxalpha = 20;

guidata(hObject, handles);

% Now compute polars for pareto front, so that smooth navigation follows
 % handles.GM.CFG.evaluate_population(handles.x_sorted);
% Make shortened version for debugging purposes...
handles.GM.CFG.evaluate_population(handles.x_sorted(2:end,:));

handles.experiment_matrix_cell = handles.GM.CFG.last_experiment_set;

% Update handles object to keep track
guidata(hObject, handles);

% LOAD CONTEXT
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Choose the case file to run
[filename, pathname, ~] = uigetfile('*.m', 'Choose case M-file');
% Add its folder to the matlab path
addpath(pathname);
% Store case name for reference
handles.casename = filename(1:end-2);

% Execute Preparation of file
eval(handles.casename);

% Define GUI handle in system context object
SC.hObject = hObject;

% Store optimization manager in workspace, together with most used objects
handles.GM = GM;
handles.SC = SC;
handles.SD = SD;

% Make Objects available in Base Workspace
assignin('base' , 'GM' , GM);

% Update handles object to keep track
guidata(hObject, handles);

% SAVE AIRFOILS TO FILE
% --- Executes on button press in savButton.
function savButton_Callback(hObject, eventdata, handles)
% hObject    handle to savButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Ask user for filename (supply image save filenmae as default
[fname, pathname] = uiputfile('*.air' , 'Select Image to Save');

% Remove filename extension and store in handles object
handles.fname_root = fname(1:end-4);

% Update handles object to keep track
guidata(hObject, handles);

disp('Hello')


% Now Extract data from table to prepare plots
uicell = get(handles.uitable1, 'Data');

% Find which airfoils should appear
plot_indices = find(cell2mat(uicell(:,1)));

save_list_of_airfoils_on_pareto(plot_indices , handles)

% SAVE IMAGES
% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Ask user for filename
[fname, pathname] = uiputfile('*.pdf' , 'Select Image to Save')

% Remove filename extension and store in handles object
handles.fname_root = fname(1:end-4);

% Update handles object to keep track
guidata(hObject, handles);

disp('Hello')


% Now Extract data from table to prepare plots
uicell = get(handles.uitable1, 'Data');

% Find which airfoils should appear
plot_indices = find(cell2mat(uicell(:,1)));

plot_list_of_airfoils_on_pareto(plot_indices , handles, true)

function  [fval_sorted , x_sorted ] = sort_pareto_data(fval , x)
    % Sort for 2nd objective
    n_sort_objective = 2;

    % Make sorting index
    [fval_sorted fval_index] = sort(fval);


    %%% Now sort phenotypes
    % Assign memory for sorted objectives array
    fval_sorted = zeros(size(fval));
    % (this can be put into a loop to handle arbitrary number of
    % objectives)
    fval_sorted(:,1) = fval(fval_index(:,n_sort_objective),1);
    fval_sorted(:,2) = fval(fval_index(:,n_sort_objective),2);

    %%% Now sort genotypes
    x_size = size(x);
    % Prealocate memory for sorted array
    x_sorted = zeros(x_size);
    for n_column = 1:x_size(2)
        x_sorted(:,n_column) = x(fval_index(:,n_sort_objective),n_column);
    end

% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

% 
disp('Hello')

% Extract data from table
uicell = get(handles.uitable1, 'Data');

% Find which airfoils should appear
plot_indices = find(cell2mat(uicell(:,1)));

% Plot them
plot_list_of_airfoils_on_pareto(plot_indices , handles, false);

function plot_list_of_airfoils_on_pareto(plot_indices , handles, for_print)

if for_print == false
    axes_airfoil = handles.axes_airfoil;
    axes_clalpha = handles.axes_clalpha;
    axes_clcd    = handles.axes_clcd;
end

if for_print == true
    f11 = figure(11);
    axes_airfoil = axes();
    f12 = figure(12);
    axes_clalpha = axes();
    f13 = figure(13);
    axes_clcd    = axes();
end

    

disp('Hello')  
% Disable Hold to restart plots
hold(axes_airfoil, 'off')
hold(axes_clalpha , 'off')
hold(axes_clcd , 'off')

% Make color code
%rainbow_list = jet(length(plot_indices));
rainbow_list = cool(length(plot_indices));


legend_list = cell(length(plot_indices),1);

for n_plot = 1:length(plot_indices) % replace by for loop when debug complete
    %Extract Current Result
    result = handles.experiment_matrix_cell{plot_indices(n_plot)};
    plot(axes_airfoil, result.coordinates.tx , result.coordinates.tz, 'Color' , rainbow_list(n_plot, :));
    
    % Now plot Clalpha
    plot(axes_clalpha , result.ap.raw_data.alpha , result.ap.raw_data.cl , 'Color' , rainbow_list(n_plot, :));
    
    % Now plot CdCl
    plot(axes_clcd, result.ap.raw_data.cd , result.ap.raw_data.cl , 'Color' , rainbow_list(n_plot, :));
       
    % Restablish holding to overlay next selecteed airfoil
    hold(axes_airfoil, 'on');
    hold(axes_clalpha, 'on');
    hold(axes_clcd   , 'on');
    
    cl_over_cd_max = max(result.ap.raw_data.cl ./ result.ap.raw_data.cd);
    thickness      = max(result.coordinates.tz) - min(result.coordinates.tz);
    
    legend_list{n_plot} = [ 'Thick = ' , num2str(thickness,3) , ' L/D_{max} =  ' , num2str(cl_over_cd_max,3)];
end


if for_print == false
% Plot upper and lower bounds of airfoil shapes
% These boudns get confusing in printing
[tx_lb_ext , tz_lb_ext] = handles.SD.generate_coordinates(200, handles.GM.CM.lb_ext);
[tx_ub_ext , tz_ub_ext] = handles.SD.generate_coordinates(200, handles.GM.CM.ub_ext);
plot(axes_airfoil , tx_lb_ext , tz_lb_ext, 'Color', [0.7 0.7 0.7])
plot(axes_airfoil , tx_ub_ext , tz_ub_ext , 'Color', [0.7 0.7 0.7]);
end

% Now polishing elements
grid(axes_airfoil   , 'on');
axis(axes_airfoil   , [-0.05 1.05 -0.25 0.25]);
xlabel(axes_airfoil , 'x/c');

axis(axes_clcd    , [0 0.03 -0.5 2]);
grid(axes_clcd    , 'on');
xlabel(axes_clcd  , 'C_{d}');
    
axis(axes_clalpha , [-10 20 -0.5 2]);
grid(axes_clalpha , 'on');
xlabel(axes_clalpha , '\alpha (deg)');

legend(axes_clcd , legend_list , 'Location' , 'Best');

if for_print == true
   % Prepare printing and save if enabled

   % Set orthonormal axes for airfoils
   axis(axes_airfoil   , [-0.05 1.05 -0.36 0.36]);

   set(f11, 'PaperType', 'A5');
   set(f12, 'PaperType', 'A5');
   set(f13, 'PaperType', 'A5');

   orient(f11, 'landscape');
   orient(f12, 'landscape');
   orient(f13, 'landscape');
   
   print(f11 , '-dpdf' , [handles.fname_root '_airfoils.pdf']);
   print(f12 , '-dpdf' , [handles.fname_root '_cl_alpha.pdf']);
   print(f13 , '-dpdf' , [handles.fname_root '_cl_cd.pdf']);
end

function save_list_of_airfoils_on_pareto(plot_indices , handles)
  
disp('Hello, I am saving')  


for n_plot = 1:length(plot_indices) % replace by for loop when debug complete
    %Extract Current Result
    result = handles.experiment_matrix_cell{plot_indices(n_plot)};
    
    % Get coordinates and compose them into an array with the right shape
    % for saving
    coordinate_array = [result.coordinates.tx , result.coordinates.tz];     %#ok<NASGU>
    % Compute thickness to write it in filename
    thickness      = max(result.coordinates.tz) - min(result.coordinates.tz);
    
    % Save !
    save([handles.fname_root,'.n_',num2str(n_plot),'.Thick', num2str(thickness,3) '.air'] , 'coordinate_array', '-ascii');
    x = result.x; %#ok<NASGU>
    save([handles.fname_root,'.n_',num2str(n_plot),'.Thick', num2str(thickness,3) '.coord'] , 'x' , '-ascii');      
end

% Done !

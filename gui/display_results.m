function display_results(hObject , results)
% hObject    handle to any object in figure
handles = guidata(hObject);


% Generate Airfoil Coorinates
[tx , tz] = handles.SD.generate_coordinates(200, results.x);
[tx_lb_ext , tz_lb_ext] = handles.SD.generate_coordinates(200, handles.GM.CM.lb_ext);
[tx_ub_ext , tz_ub_ext] = handles.SD.generate_coordinates(200, handles.GM.CM.ub_ext);


% Plot Airfoil together with upper and lower bounds
plot(handles.axes_airfoil   , tx, tz)
hold(handles.axes_airfoil , 'on')
plot(handles.axes_airfoil , tx_lb_ext , tz_lb_ext, 'Color', [0.7 0.7 0.7])
plot(handles.axes_airfoil , tx_ub_ext , tz_ub_ext , 'Color', [0.7 0.7 0.7]);
hold(handles.axes_airfoil , 'off')
grid(handles.axes_airfoil   , 'on');
axis(handles.axes_airfoil   , [-0.05 1.05 -0.25 0.25]);
xlabel(handles.axes_airfoil , 'x/c');

% Now!
% Extract Polar data
if strcmp(class(results.ap) , 'aerodynamic_polar')
    % If data is valid
    alpha = results.ap.raw_data.alpha;
    cl    = results.ap.raw_data.cl;
    cd    = results.ap.raw_data.cd;
    
    % Write description of what's going on
    % Proceed to write description of airfoil properties
    [ld_max i_ld_max]   = max(cl ./ cd);
    cl_ld_max           = cl(i_ld_max);
    alpha_stall         = results.ap.alpha_cl_max_local_limited;
    cl_stall            = results.ap.cl_max_local_limited;
    
    % Get Information on Cost function Values (only if we are on single_sim
    % case!)
    if length(handles.GM.CFG.simulation_worker_list) == 1
        % Single Sim case can be interpreted
        val_obj         = handles.GM.CFG.interpret_sample({results});
    else
        % Multi Sim case cannot be interpreted with this approach! Return
        % to this later! With a clean code this time, would be nice!!!
        val_obj         = [0 0];
    end
    
    data_line       = [ld_max cl_ld_max cl_stall alpha_stall val_obj(1) val_obj(2)];
    
    
    % Obtain former data
    table_data_cell = get(handles.current_airfoil_table , 'Data');
    if isempty(table_data_cell{1,1})
        % If we are in initialization run define max and min
        max_line = data_line;
        min_line = data_line;
    else
        % Otherwise, update max and min
        table_data_array = cell2mat(table_data_cell);
        % Obtain current extreme values
        max_line = table_data_array(2,:);
        min_line = table_data_array(3,:);
        % Compare and update if applicable
        max_line = max([ data_line ; max_line]);
        min_line = min([ data_line ; min_line]);
    end
    
    % Recompose table
    table_data_array    = [data_line ; max_line ; min_line];
    table_data_cell = mat2cell(table_data_array,ones(1 , size(table_data_array) * [1;0]),ones(1 , size(table_data_array) * [0;1])); %#ok<MMTC>
    % Draw!
    set(handles.current_airfoil_table, 'Data' , table_data_cell);
    
    % Plot Polar
    
    % Determine range over which to plot  (Keep track of maximum cl over which
    % so plot does not change scale with each profile, but increases scale when
    % necessary)
    if isfield(handles , 'maxcl')
        handles.maxcl       = max(max(cl)   , handles.maxcl);
        handles.mincl       = min(min(cl)   , handles.mincl);
        
        handles.maxalpha    = max(max(alpha), handles.maxalpha);
        handles.minalpha    = min(min(alpha), handles.minalpha);
    else
        % In case adaptive plot bounds were not initialized, start things up!
        handles.maxcl       = max(max(cl), 0.8);
        handles.mincl       = min(min(cl), -0.2);
        
        handles.maxalpha    = max(max(alpha), 10);
        handles.minalpha    = min(min(alpha), 0);
    end
    
    % Plot and set axis!
    plot(handles.axes_clcd    , cd    , cl);
    axis(handles.axes_clcd    , [0 0.2 handles.mincl 1.1*handles.maxcl]);
    grid(handles.axes_clcd    , 'on');
    xlabel(handles.axes_clcd  , 'C_{l}');
    
    plot(handles.axes_clalpha , alpha , cl);
    axis(handles.axes_clalpha , [handles.minalpha handles.maxalpha handles.mincl 1.1*handles.maxcl]);
    grid(handles.axes_clalpha , 'on');
    xlabel(handles.axes_clalpha , '\alpha (deg)');
end

if ~strcmp(class(results.ap) , 'aerodynamic_polar')
    try
    % Polar data is invalid (write and plot 0's)
    plot(handles.axes_clcd    , 0    , 0);
    axis(handles.axes_clcd    , [0 0.2 handles.mincl 1.1*handles.maxcl]);
    grid(handles.axes_clcd    , 'on');
    xlabel(handles.axes_clcd  , 'C_{l}');
    
    plot(handles.axes_clalpha , 0 , 0);
    axis(handles.axes_clalpha , [handles.minalpha handles.maxalpha handles.mincl 1.1*handles.maxcl]);
    grid(handles.axes_clalpha , 'on');
    xlabel(handles.axes_clalpha , '\alpha (deg)');
    
    
    table_data_cell = get(handles.current_airfoil_table , 'Data');
    % Get min and max, and modify data of table to show 0's at current data
    table_data_array = cell2mat(table_data_cell);
    % Obtain current extreme values
    max_line = table_data_array(2,:);
    min_line = table_data_array(3,:);    
    data_line = zeros(size(max_line));
    table_data_array    = [data_line ; max_line ; min_line];
    table_data_cell = mat2cell(table_data_array,ones(1 , size(table_data_array) * [1;0]),ones(1 , size(table_data_array) * [0;1])); %#ok<MMTC>
    % Draw!
    set(handles.current_airfoil_table, 'Data' , table_data_cell);
    catch
        disp('Failed Fallback in display_results')
    end
end


% Update handles object to keep track
guidata(hObject, handles);


% BL data reader

% Set file to read
filename_root       = 'du25plr16';
extension_cfx_file   = '.cfx';
extension_cpx_file   = '.cpx';
extension_dst_file   = '.dst';
extension_tet_file   = '.tet';

% Load file data (this is all numbers!)
raw_bldata_cfx = load([filename_root extension_cfx_file]);
raw_bldata_cpx = load([filename_root extension_cpx_file]);
raw_bldata_dst = load([filename_root extension_dst_file]);
raw_bldata_tet = load([filename_root extension_tet_file]);

% Extract data blocks from raw_bldata arrays
base_bldata_aoa  = raw_bldata_cfx(1    ,2:end);     % Angle of attack data

base_bldata_cfx  = raw_bldata_cfx(3:end,2:end);     % Skin friction data
base_bldata_cpx  = raw_bldata_cpx(3:end,2:end);     % Pressure coefficient data
base_bldata_dst  = raw_bldata_dst(3:end,2:end);     % Displacement thickness data
base_bldata_tet  = raw_bldata_tet(3:end,2:end);     % Momentum thickness (over chord) data 

base_bldata_xun  = raw_bldata_cfx(3:end,1    );     % x position data

% Sort AoA Index (for bldata)
[sorted_aoa , sort_index_aoa] = sort(base_bldata_aoa);

% Identify and index unique AoA values (for bldata)
[unique_aoa , unique_index_aoa, ic_bldata_aoa ] = unique(sorted_aoa);

% Index for unique and sorted AoA values (for bldata)
unique_sorting_index_aoa = sort_index_aoa(unique_index_aoa);

% Regroup and sort free data fields
unique_sorted_bldata_cfx = base_bldata_cfx(:,unique_sorting_index_aoa);
unique_sorted_bldata_cpx = base_bldata_cpx(:,unique_sorting_index_aoa);
unique_sorted_bldata_dst = base_bldata_dst(:,unique_sorting_index_aoa);
unique_sorted_bldata_tet = base_bldata_tet(:,unique_sorting_index_aoa);

% Now create a monotonically growing x coordinate
monotonic_xun = [0 ; cumsum(abs(diff(base_bldata_xun)))];
% And make sure it is scaled to avoid numerical artifacts
monotonic_xun = monotonic_xun / monotonic_xun(end) * 2;

% Create interpolation grids (independent variables)
[aoa_grid, xun_grid] = meshgrid(unique_aoa, monotonic_xun);

% Create base interpolant functions
base_cfx_fun = @(aoa, xun) interp2(aoa_grid, xun_grid, unique_sorted_bldata_cfx, aoa, xun);
base_cpx_fun = @(aoa, xun) interp2(aoa_grid, xun_grid, unique_sorted_bldata_cpx, aoa, xun);
base_dst_fun = @(aoa, xun) interp2(aoa_grid, xun_grid, unique_sorted_bldata_dst, aoa, xun);
base_tet_fun = @(aoa, xun) interp2(aoa_grid, xun_grid, unique_sorted_bldata_tet, aoa, xun);


% Some plotting (Cp)
aoa = 8;
xun_plot = linspace(1, 0);
plot(xun_plot, -base_cpx_fun(aoa, 1-xun_plot));
hold on; grid on;
plot(  xun_plot, -base_cpx_fun(aoa, xun_plot+1));
legend('Upper Side', 'Lower Side');
% Some plotting (Cf)
aoa = 8;
plot(xun_plot, base_cfx_fun(aoa, 1-xun_plot));
hold on; grid on;
plot(  xun_plot, base_cfx_fun(aoa, xun_plot+1));
legend('Upper Side', 'Lower Side');

% Back to work
% And sided interpolant functions
cpx_fun_top  = @(aoa, xun) -base_cpx_fun(aoa, 1-xun  );
cpx_fun_bot  = @(aoa, xun) -base_cpx_fun(aoa,   xun+1);

cfx_fun_top  = @(aoa, xun)  base_cfx_fun(aoa, 1-xun  );
cfx_fun_bot  = @(aoa, xun)  base_cfx_fun(aoa,   xun+1);

cfx_dst_top  = @(aoa, xun)  base_dst_fun(aoa, 1-xun  );
cfx_dst_bot  = @(aoa, xun)  base_dst_fun(aoa,   xun+1);

cfx_tet_top  = @(aoa, xun)  base_tet_fun(aoa, 1-xun  );
cfx_tet_bot  = @(aoa, xun)  base_tet_fun(aoa,   xun+1);

% New Plotting!
aoa = 8;
xun_plot = linspace(1, 0, 1000);
plot(xun_plot, cpx_fun_top(aoa, xun_plot));
hold on; grid on;
plot(  xun_plot, cpx_fun_bot(aoa, xun_plot));
legend('Upper Side', 'Lower Side');
% Validation
plot(base_bldata_xun, -unique_sorted_bldata_cpx( :, 61), 'x')



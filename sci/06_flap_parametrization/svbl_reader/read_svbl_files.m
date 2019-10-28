function [bldata] = read_svbl_files(filename_root)
%READ_SVBL_FILES - BL data reader for Rfoil SVBL files
%   Reads data from four boundary layer data files (cfx, cpx, dst and tet
%   extensions).  Based on bl_reader_alfa3

% Define bl_data as structure
bldata               = struct();

% Define Extensions of SVBL files
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
base_bldata_isp  = raw_bldata_cfx(2    ,2:end);     % Index of stagnation point

base_bldata_cfx  = raw_bldata_cfx(3:end,2:end);     % Skin friction data
base_bldata_cpx  = raw_bldata_cpx(3:end,2:end);     % Pressure coefficient data
base_bldata_dst  = raw_bldata_dst(3:end,2:end);     % Displacement thickness data
base_bldata_tet  = raw_bldata_tet(3:end,2:end);     % Momentum thickness (over chord) data 

base_bldata_xun  = raw_bldata_cfx(3:end,1    );     % x position data

% Sort AoA Index (for bldata)
[sorted_aoa , sort_index_aoa] = sort(base_bldata_aoa);

% Identify and index unique AoA values (for bldata)
[unique_sorted_aoa , unique_index_aoa, ~ ] = unique(sorted_aoa);

% Index for unique and sorted AoA values (for bldata)
unique_sorting_index_aoa = sort_index_aoa(unique_index_aoa);

% Regroup and sort free data fields
unique_sorted_bldata_isp = base_bldata_isp(:,unique_sorting_index_aoa);
unique_sorted_bldata_cfx = base_bldata_cfx(:,unique_sorting_index_aoa);
unique_sorted_bldata_cpx = base_bldata_cpx(:,unique_sorting_index_aoa);
unique_sorted_bldata_dst = base_bldata_dst(:,unique_sorting_index_aoa);
unique_sorted_bldata_tet = base_bldata_tet(:,unique_sorting_index_aoa);

% Now create a monotonically growing x coordinate
monotonic_xun = [0 ; cumsum(abs(diff(base_bldata_xun)))];
% And make sure it is scaled to avoid numerical artifacts
monotonic_xun = monotonic_xun / monotonic_xun(end) * 2;

% Create interpolation grids (independent variables)
[aoa_grid, xun_grid] = meshgrid(unique_sorted_aoa, monotonic_xun);

% Create base interpolant functions
base_cfx_fun = @(aoa, xun) interp2(aoa_grid, xun_grid, unique_sorted_bldata_cfx, aoa, xun);
base_cpx_fun = @(aoa, xun) interp2(aoa_grid, xun_grid, unique_sorted_bldata_cpx, aoa, xun);
base_dst_fun = @(aoa, xun) interp2(aoa_grid, xun_grid, unique_sorted_bldata_dst, aoa, xun);
base_tet_fun = @(aoa, xun) interp2(aoa_grid, xun_grid, unique_sorted_bldata_tet, aoa, xun);


% Now, document position of stagnation point
% Find position of leading edge
[~, ile] = min(base_bldata_xun);
% Make array for side of stagnation point as function of angle of attack
% (rfoil convention: 1 = upper side, 2 = lower side)
unique_sorted_bldata_isp_side = 2 - (unique_sorted_bldata_isp < ile);
% Make array for position of stagnation point as function of angle of
% attack (not smooth! means finer panelling would be useful for high AOA! check it with a plot!)
unique_sorted_bldata_xsp = transpose(base_bldata_xun(unique_sorted_bldata_isp));
% Then make interpolants for position, index and side of stagnation point
bldata.base_xsp_fun      = @(aoa) interp1(unique_sorted_aoa,unique_sorted_bldata_xsp, aoa);                  % For position of stagnation point
bldata.base_isp_fun      = @(aoa) interp1(unique_sorted_aoa,unique_sorted_bldata_isp, aoa);                  % For index of stagnation point
bldata.base_isp_side_fun = @(aoa) interp1(unique_sorted_aoa,unique_sorted_bldata_isp_side, aoa, 'nearest');  % For side of stagnation point (nearest, to avoid non integer values!)

% And sided interpolant functions for boundary layer quantities
bldata.cpx_fun_top  = @(aoa, xun) -base_cpx_fun(aoa, 1-xun);
bldata.cpx_fun_bot  = @(aoa, xun) -base_cpx_fun(aoa, 1+xun);

bldata.cfx_fun_top  = @(aoa, xun)  base_cfx_fun(aoa, 1-xun);
bldata.cfx_fun_bot  = @(aoa, xun)  base_cfx_fun(aoa, 1+xun);

bldata.dst_fun_top  = @(aoa, xun)  base_dst_fun(aoa, 1-xun);
bldata.dst_fun_bot  = @(aoa, xun)  base_dst_fun(aoa, 1+xun);

bldata.tet_fun_top  = @(aoa, xun)  base_tet_fun(aoa, 1-xun);
bldata.tet_fun_bot  = @(aoa, xun)  base_tet_fun(aoa, 1+xun);

% Add a sided interpolant function for shape factors
bldata.h_fun_top    = @(aoa, xun)  base_dst_fun(aoa, 1-xun) ./ base_tet_fun(aoa, 1-xun);
bldata.h_fun_bot    = @(aoa, xun)  base_dst_fun(aoa, 1+xun) ./ base_tet_fun(aoa, 1+xun);

% Also add angle of attack and monotonic_xun ranges
bldata.aoa_range   = unique_sorted_aoa;
bldata.xun_top_range = base_bldata_xun(ile:(-1):1  );
bldata.xun_bot_range = base_bldata_xun(ile:(+1):end);

% And angle of attack and monotonic_xun grids
bldata.aoa_grid    = aoa_grid;
bldata.xun_grid    = xun_grid;
% Finally, add raw xun stes (monotonic and base)
bldata.xun_range   = monotonic_xun;
bldata.base_xun    = base_bldata_xun;


end


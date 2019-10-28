% Script for Reading a Polar File from TR824
% 
%   alfa5 tested manually on 131.ALL and 132.ALL and exhaustively for
%   chapters A and B. Chapter has slightly different filename conventions
%   but we will not use that data explicitly.


% % Now read some data into the system

% Set Name of File to Read
filename = './data/experimental/NACA/TR824-Digitized/132.ALL';

[datasets] = TR824_reader.read_datasets_from_file(filename);

% Let us now extract a tripplet of datasets
alcl_polar_dataset     = datasets.A01_data;
alcm_polar_dataset     = datasets.A07_data;
clcd_polar_dataset     = datasets.B01_data;

% % Take care of Cl polar
% Extract data from cl_polar_dataset
al_alcl_polar_raw      = alcl_polar_dataset(:,1);
cl_alcl_polar_raw      = alcl_polar_dataset(:,2);

% Sort cl_polar dataset
[~, sort_index_cl_polar] =  sort(al_alcl_polar_raw);
al_alcl_polar_sorted   = al_alcl_polar_raw(sort_index_cl_polar);
cl_alcl_polar_sorted   = cl_alcl_polar_raw(sort_index_cl_polar);

% Now find monotonous streak
% Compute forward difference of Cl vs alpha
dcl_alcl_polar_sorted  = diff(cl_alcl_polar_sorted);
% Find indices of points where polar is decreasing
index_dcl_decreasing   = find(dcl_alcl_polar_sorted < 0);
% Find a point near zero angle of attack
[~, index_al_min_alcl_polar_sorted] = min(al_alcl_polar_sorted.^2);
% Find highest aoa below near-zero piont where polar is decreasing 
index_start_monotonous = max(index_dcl_decreasing(index_dcl_decreasing < index_al_min_alcl_polar_sorted)) + 1;
% Find lowest aoa above near-zero piont where polar is decreasing 
index_end_monotonous   = min(index_dcl_decreasing(index_dcl_decreasing > index_al_min_alcl_polar_sorted));

% Now make monotonous streaks of cl polar
al_alcl_polar_monotonous = al_alcl_polar_sorted(index_start_monotonous:index_end_monotonous);
cl_alcl_polar_monotonous = cl_alcl_polar_sorted(index_start_monotonous:index_end_monotonous);

% Now interpolants on monotonous streaks of cl polar
al_alcl_polar_fun      = @(cl) interp1(cl_alcl_polar_monotonous , al_alcl_polar_monotonous, cl );

% % Take care of Cd polar
% Now extract data from cd_polar_dataset (cl, cd)
cl_clcd_polar_raw      = clcd_polar_dataset(:,1);
cd_clcd_polar_raw      = clcd_polar_dataset(:,2);

% Now sort cl, cd polar (not used for reinterpolation)
[~, sort_index_clcd_polar] =  sort(cl_clcd_polar_raw);
cl_clcd_polar_sorted   = cl_clcd_polar_raw(sort_index_clcd_polar);
cd_clcd_polar_sorted   = cd_clcd_polar_raw(sort_index_clcd_polar);


% Now reinterpolate alpha from monotonous streak of cl polar to find
% alpha's of cd polar
al_reinterpolated      = al_alcl_polar_fun(cl_clcd_polar_raw);
% Filter NaNs out
al_filtered            = al_reinterpolated(not(isnan(al_reinterpolated)));
cl_filtered            = cl_clcd_polar_raw(not(isnan(al_reinterpolated)));
cd_filtered            = cd_clcd_polar_raw(not(isnan(al_reinterpolated)));
% Sort 
[~, index_al_filtered] = sort(al_filtered);

% Make alpha, cl, cd polar data
al_alclcd_polar         = al_filtered(index_al_filtered);
cl_alclcd_polar         = cl_filtered(index_al_filtered);
cd_alclcd_polar         = cd_filtered(index_al_filtered);

% % Finally, take care of Cm polar
% Extract data from cm_polar_dataset
al_alcm_polar_raw      = alcm_polar_dataset(:,1);
cm_alcm_polar_raw      = alcm_polar_dataset(:,2);

% Sort cl_polar dataset
[~, sort_index_cl_polar] =  sort(al_alcm_polar_raw);
al_alcm_polar_sorted   = al_alcm_polar_raw(sort_index_cl_polar);
cm_alcm_polar_sorted   = cm_alcm_polar_raw(sort_index_cl_polar);

% Now store data into a structure
% First : alpha, cl polar
processed_polars.al_alcl_polar    = al_alcl_polar_sorted;
processed_polars.cl_alcl_polar    = cl_alcl_polar_sorted;
% Second: cl, cd polar
processed_polars.cl_clcd_polar    = cl_clcd_polar_sorted;
processed_polars.cd_clcd_polar    = cd_clcd_polar_sorted;
% Third : alpha, cl, cd polar
processed_polars.al_alclcd_polar  = al_alclcd_polar;
processed_polars.cl_alclcd_polar  = cl_alclcd_polar;
processed_polars.cd_alclcd_polar  = cd_alclcd_polar;
% Fourth: alpha, cm polar
processed_polars.al_alcm_polar    = al_alcm_polar_sorted;
processed_polars.cm_alcm_polar    = cm_alcm_polar_sorted;






% Plot diagnostics
figure(1)
plot(al_alcl_polar_sorted    , cl_alcl_polar_sorted    , 'x-'); hold on;
plot(al_alcl_polar_monotonous, cl_alcl_polar_monotonous, 'o-'); grid on;
plot(al_filtered             , cl_filtered             , '*-');

figure(2)
subplot(122)
plot(processed_polars.al_alcl_polar   , processed_polars.cl_alcl_polar   , 'x-'); hold on;
plot(processed_polars.al_alclcd_polar , processed_polars.cl_alclcd_polar , 'o-'); grid on;
xlabel('\alpha (deg)'); ylabel('C_l');
subplot(121)
plot(processed_polars.cd_clcd_polar   , processed_polars.cl_clcd_polar   , 'x-'); hold on;
plot(processed_polars.cd_alclcd_polar , processed_polars.cl_alclcd_polar , 'o-'); grid on;
xlabel('C_d')         ; ylabel('C_l');
legend('Sorted Data' , 'Reinterpolated Data', 'Location', 'East')


% Create interpolants 







% 
%Re_dataset_01 = 3e6;
%cl_polar_alpha_data_dataset_01

%dataset_01_alpha_cl = ;
%dataset_01_alpha_for_cl = ;
%dataset_01_alpha_for_cl = ;


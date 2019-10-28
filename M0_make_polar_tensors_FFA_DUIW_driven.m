% Merges FFA and DUWP into polar tensor, using externally provided fID
% e.g. fID = 'S00E35';

%% Free case
% Load FFA airfoil polars (for free case)
DUWP210_free  = load(['rotor_integration/airfoil_families/DU-IW/' , fID , '_polar_mesh_free.mat']);
FFAw3241_free = load('rotor_integration/airfoil_families/FFA/FFAw3241_polar_mesh_free.mat');
FFAw3301_free = load('rotor_integration/airfoil_families/FFA/FFAw3301_polar_mesh_free.mat');
FFAw3360_free = load('rotor_integration/airfoil_families/FFA/FFAw3360_polar_mesh_free.mat');
FFAw3480_free = load('rotor_integration/airfoil_families/FFA/FFAw3480_polar_mesh_free.mat');
FFAw3600_free = load('rotor_integration/airfoil_families/FFA/FFAw3600_polar_mesh_free.mat');
%cylinder_free = load('rotor_integration/airfoil_families/FFA/cylinder_polar_mesh_free.mat');


% Make relative thickness range
tc_range    = [0.210, 0.241, 0.301, 0.360, 0.480, 0.600];
% And corresponding cell array
polar_mesh_cell = cell(4,1);
polar_mesh_cell{1} = DUWP210_free.polar_mesh;
polar_mesh_cell{2} = FFAw3241_free.polar_mesh;
polar_mesh_cell{3} = FFAw3301_free.polar_mesh;
polar_mesh_cell{4} = FFAw3360_free.polar_mesh;
polar_mesh_cell{5} = FFAw3480_free.polar_mesh;
polar_mesh_cell{6} = FFAw3600_free.polar_mesh;
% polar_mesh_cell{5} = cylinder_free.polar_mesh;
% Make polar tensors
polar_tensors = make_polar_tensors(tc_range, polar_mesh_cell); %#ok<NASGU>
% Save to file
save(['rotor_integration/airfoil_families/DU-IW/' , fID , '_free.mat'], 'polar_tensors')

%% Tripped case
% Load FFA airfoil polars (for free case)
DUWP210_trip = load(['rotor_integration/airfoil_families/DU-IW/', fID , '_polar_mesh_trip.mat']);
FFAw3241_trip = load('rotor_integration/airfoil_families/FFA/FFAw3241_polar_mesh_trip.mat');
FFAw3301_trip = load('rotor_integration/airfoil_families/FFA/FFAw3301_polar_mesh_trip.mat');
FFAw3360_trip = load('rotor_integration/airfoil_families/FFA/FFAw3360_polar_mesh_trip.mat');
FFAw3480_trip = load('rotor_integration/airfoil_families/FFA/FFAw3480_polar_mesh_trip.mat');
FFAw3600_trip = load('rotor_integration/airfoil_families/FFA/FFAw3600_polar_mesh_trip.mat');
%cylinder_trip = load('rotor_integration/airfoil_families/FFA/cylinder_polar_mesh_trip.mat');

% Make relative thickness range
tc_range    = [0.210, 0.241, 0.301, 0.360, 0.480, 0.600];
% And corresponding cell array
polar_mesh_cell = cell(4,1);
polar_mesh_cell{1} = DUWP210_trip.polar_mesh;
polar_mesh_cell{2} = FFAw3241_trip.polar_mesh;
polar_mesh_cell{3} = FFAw3301_trip.polar_mesh;
polar_mesh_cell{4} = FFAw3360_trip.polar_mesh;
polar_mesh_cell{5} = FFAw3480_trip.polar_mesh;
polar_mesh_cell{6} = FFAw3600_trip.polar_mesh;
% polar_mesh_cell{5} = cylinder_trip.polar_mesh;
% Make polar tensors
polar_tensors = make_polar_tensors(tc_range, polar_mesh_cell);
% Save to file
save(['rotor_integration/airfoil_families/DU-IW/' , fID , '_trip.mat'], 'polar_tensors')


%% Function to make polar tensors
function polar_tensors = make_polar_tensors(tc_range, polar_mesh_cell)
    % Extract angle of attack and Re ranges
    alpha_range = polar_mesh_cell{1}.alpha_range;
    Re_range    = polar_mesh_cell{2}.Re_range   ;
    
    % Now make interpolation tensor
    [alpha_tensor, Re_tensor, tc_tensor] = ndgrid(alpha_range, Re_range, tc_range);
    
    % Allocate coefficient tensors
    cl_tensor = zeros(size(alpha_tensor));
    cd_tensor = zeros(size(alpha_tensor));
    %cm_tensor = zeros(size(alpha_tensor));
    
    % Fill them in
    for n_tc = 1:length(tc_range)
        cl_tensor(:,:,n_tc) = polar_mesh_cell{n_tc}.cl_mesh;
        cd_tensor(:,:,n_tc) = polar_mesh_cell{n_tc}.cd_mesh;
        %cm_tensor(:,:,n_tc) = polar_mesh_cell{n_tc}.cm_mesh;
    end
    
%     % And make interpolating function
%     cl_tensor_fun = @(alpha, Re, tc) interpn(alpha_tensor, Re_tensor, tc_tensor, cl_tensor,  alpha, Re, tc);
%     cd_tensor_fun = @(alpha, Re, tc) interpn(alpha_tensor, Re_tensor, tc_tensor, cl_tensor,  alpha, Re, tc);
%     cm_tensor_fun = @(alpha, Re, tc) interpn(alpha_tensor, Re_tensor, tc_tensor, cl_tensor,  alpha, Re, tc);
    
    % Bundle all into a structure
    polar_tensors = struct();
    % Independent variable ranges
    polar_tensors.alpha_range   = alpha_range;
    polar_tensors.Re_range      = Re_range;
    polar_tensors.tc_range      = tc_range;
    % Independent variable tensors
    polar_tensors.alpha_tensor  = alpha_tensor;
    polar_tensors.Re_tensor     = Re_tensor;
    polar_tensors.tc_tensor     = tc_tensor;
    % Dependent variable tensors
    polar_tensors.cl_tensor     = cl_tensor;
    polar_tensors.cd_tensor     = cd_tensor;
    %polar_tensors.cm_tensor     = cm_tensor;
%     % Function handles
%     polar_tensors.cl_tensor_fun = cl_tensor_fun;
%     polar_tensors.cd_tensor_fun = cd_tensor_fun;
%     polar_tensors.cm_tensor_fun = cm_tensor_fun;
end





function data_struct = read_dat_file_into_data_struct(filename)

% Read File Header
header = read_header_from_dat_file(filename);

% Process Header to get variable names
var_names = header(2,:);

% Process to estimate array size
I = str2num(header{3,2});
J = str2num(header{3,4});

% Process to estimate array size
startRow = 4;
endRow   = startRow + I*J -1; % Last row is only zeros!

% Read data array
data_array = read_array_from_dat_file(filename, startRow, endRow);

% Put Data into a Structure
data_struct = struct();

data_struct.z              =   reshape(data_array(:,1), I , J)' * 1e-3; % Get it into meters
data_struct.y              =   reshape(data_array(:,2), I , J)' * 1e-3; % Get it into meters

data_struct.w              =   reshape(data_array(:,3), I , J)';       % 
data_struct.v              =   reshape(data_array(:,4), I , J)';       %
data_struct.u              = - reshape(data_array(:,5), I , J)';       % Make free-stream positive!
data_struct.uvw_length     =   reshape(data_array(:,6), I , J)';       %

data_struct.w_rms          =   reshape(data_array(:,7), I , J)';
data_struct.v_rms          =   reshape(data_array(:,8), I , J)';
data_struct.u_rms          =   reshape(data_array(:,9), I , J)';
data_struct.uvw_rms_length =   reshape(data_array(:,10), I , J)';

data_struct.rst_zy         =   reshape(data_array(:,11), I , J)';
data_struct.rst_zx         =   reshape(data_array(:,12), I , J)';
data_struct.rst_yx         =   reshape(data_array(:,13), I , J)';

data_struct.rst_zz         =   reshape(data_array(:,14), I , J)';
data_struct.rst_yy         =   reshape(data_array(:,15), I , J)';
data_struct.rst_xx         =   reshape(data_array(:,16), I , J)';
data_struct.I              =   J;
data_struct.J              =   I;

% Make perturbation field estimate
u_pert = zeros(size(data_struct.u));
for i_line = 1:data_struct.I
    u_pert(i_line,:) = data_struct.u(i_line,:) - mean(data_struct.u(i_line,:));
end
data_struct.u_pert         = u_pert;

% Make interpolant functions
data_struct.u_fun          = @(z,y) interp2(data_struct.z, data_struct.y, data_struct.u, z , y);
data_struct.v_fun          = @(z,y) interp2(data_struct.z, data_struct.y, data_struct.v, z , y);
data_struct.w_fun          = @(z,y) interp2(data_struct.z, data_struct.y, data_struct.w, z , y);

data_struct.interp_fun     = @(field,z,y) interp2(data_struct.z, data_struct.y, field, z , y);

end






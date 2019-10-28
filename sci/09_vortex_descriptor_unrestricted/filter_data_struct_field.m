function data_struct = filter_data_struct_field( data_struct, fieldname, z_min, z_max, y_min, y_max)
% Filter out unreliable data points in data_struct, average in spanwise
% direction

% Get info on filtering window
[~, z_min_index] = min((data_struct.z(1,:) - z_min).^2);
[~, z_max_index] = min((data_struct.z(1,:) - z_max).^2);
[~, y_min_index] = min((data_struct.y(:,1) - y_min).^2);
[~, y_max_index] = min((data_struct.y(:,1) - y_max).^2);

% Filter out data points!
uf_field = data_struct.(fieldname);
f_field  = uf_field;
for y_index = min(y_min_index,y_max_index):max(y_min_index,y_max_index)
    uf_field_z_min = uf_field(y_index, z_min_index);
    uf_field_z_max = uf_field(y_index, z_max_index);
    z_min_current =  data_struct.z(y_index, z_min_index);
    z_max_current =  data_struct.z(y_index, z_max_index);
    for z_index = min(z_min_index,z_max_index):max(z_min_index,z_max_index)
        z_current     =  data_struct.z(y_index, z_index);
        f_field(y_index, z_index) = uf_field_z_min + ...
            (uf_field_z_max-uf_field_z_min) / (z_max_current-z_min_current) ...
                                            * (z_current - z_min_current);
        %f_field(y_index, z_index) = uf_field_z_min;
        
    end
end

% Update data struct with filtered field before returning it!
data_struct.(fieldname) = f_field;

end


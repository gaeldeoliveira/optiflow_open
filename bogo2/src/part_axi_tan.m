function [var_axi, var_tan] = part_axi_tan(var_serial)

    var_axi = var_serial(1:(length(var_serial)/2));
    var_tan = var_serial((length(var_serial)/2+1):end);
end
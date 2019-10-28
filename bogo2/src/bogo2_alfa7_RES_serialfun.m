function [RES_serial] = bogo2_alfa7_RES_serialfun(BC, BD, a_serial)
% Serializes bogo2_alfa6_RES_fun inputs and outputs
    
    % Part primary variables
    [a_axi, a_tan] = part_axi_tan(a_serial);
    
    % Run function
    [RES_axi, RES_tan] = bogo2_alfa7_RES_fun(BC, BD, a_axi, a_tan);
    
    % Serialize outputs
    RES_serial = serialize_axi_tan(RES_axi, RES_tan);

end
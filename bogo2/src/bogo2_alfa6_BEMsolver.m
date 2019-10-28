function [RES_axi, RES_tan, BD, BR] = bogo2_alfa6_BEMsolver(BC)
        % % Make discretization
        [~, ~, BD, ~] = bogo2_alfa6_RES_fun(BC, [], [], []);
        % Extract initial guess
        a_serial_0 = serialize_axi_tan(BD.a_axi_0, BD.a_tan_0);
        % Make serial function wrapper
        res_serial_fun = @(a_serial) bogo2_alfa6_RES_serialfun(BC, BD, a_serial);
        
        % % Now solve with non-linear solution algorithm
        % Make options for solution algorithm
        options = optimoptions('fsolve');
        options.StepTolerance           = 1e-18;
        options.MaxFunctionEvaluations  = 1e5;
        options.FunctionTolerance       = 1e-9;
        options.OptimalityTolerance     = 1e-9;
        options.MaxIterations           = 1e5;
        options.Display                 = 'off';
        % Solve!
        a_serial                        = fsolve(res_serial_fun, a_serial_0, options);
        
        % % Documentation round
        [a_axi, a_tan] = part_axi_tan(a_serial);
        [RES_axi, RES_tan, BD, BR] = bogo2_alfa6_RES_fun(BC, BD, a_axi, a_tan);
end
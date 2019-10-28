classdef plasma_actuator_parameters
    %PLASMA_ACTUATOR_PARAMETERS Describes everything about a single plasma
    % actuator
    % 
    % This is not a handle class, but rather a simple container with
    % default parameters, to be used like a template! 
    %
    % Meant to be used as an input to the simulation protocol, with 
    %   single_plasma  - single plasma actuator, defined via
    %                     plasma_actuator_parameters structure with 
    %                     four fields:
    %                      .tpp  = plasma force field thickness over chord
    %                      .lpp  = plasma force field lenght over chord
    %                      .xpp  = plasma beginning (leading edge) over chord
    %                      .cftp = plasma force field coefficient (cftp = FTP [N] / (0.5*rho*U_inf^2*chord))
    %                      .ispp = plasma actuator side (integer, 1 = top,2 = bottom)
    %                      .cfc  = plasma skin friction correction factor
    %                      .ceip = plasma energy interaction coefficient correction factor
    %
    % Because this is not a handle subclass, it lends itself well to be used 
    % with a simulation_protocol_dynamizer in parallel computation cases 
    % (the global nature of the handle superclass would induce race conditions).
    
    
    properties
      % Actuator Parameters (some usual values are provided as default!)
      tpp  =  0.0060  % plasma force field thickness over chord
      lpp  =  0.1950  % plasma force field lenght over chord
      xpp  =  0.6500  % plasma beginning (leading edge) over chord
      cftp = -0.0012  % plasma force field coefficient (cftp = FTP [N] / (0.5*rho*U_inf^2*chord))
      ispp =  1       % plasma actuator side (integer, 1 = top,2 = bottom)
      
      % Correction factors
      cfc  = 1.0      % plasma skin friction correction factor
      ceip = 1.0      % plasma energy interaction coefficient correction factor
    end
    
    methods
        function PL = plasma_actuator_parameters()
           % Dummy Constructor function (declared explicitly!) 
        end
    end
    
end


% A simple modular Panel Code based on a free interpretation of the
% Hess-Smith method in velocity components
% Gael de Oliveira
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%       ActiHS  :   A simple modular Panel Code based on a free           %
%                   interpretation of the Hess-Smith method in velocity   %
%                   components                                            %
%                                                                         %
%       Usage   :   Standalone with script, for sail optimization, within %
%                   the kirikou-dogoro actuator codes or other codes and  %
%                   derivatives from the author                           %
%                                                                         %
%       Date    :   April 2011 to March 2017                              %
%       Author  :   Gael de Oliveira                                      %
%                                                                         %
%       License :   MIT, as the rest of this repository 		    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Load Airfoil Geometry
coord = load('airfoils/naca0012_160.air');

% Extract original coordinates
px = [coord(:,1)];
py = [coord(:,2)];

% Get Index
i_start = [1 ;length(px)+1];
i_end   = [length(px) ;2*length(px)];



% Plot
plot(px, py , 'o-')
grid on
axis equal

%% Duplicate Airfoil
ipc = inviscid_panel_case_multielement([px ; px + 2] , [py ; py], [1 ;length(px)+1], [length(px) ;2*length(px)]);

%% Create Case
ipc.rotation = false;                   % Set free stream to be straight (usual) (rotation=true introduces a correction for curved streamlines, as in VAWT applications!)
ipc.alpha = 5*pi/180;                   % Set angle of attack in radians
ipc.generate_solution;                  % Generate a solution!

% And plot results!
% Here we plot Cp but there is also some more pre-postprocessed data! 
plot(ipc.px_middle , ipc.cp_plot)
grid on
xlabel('x/c')
ylabel('Cp')
title('Inviscid Cp Plot')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Kirikou -   A simple Vorticity Equation Solver Actuator Disk Flows    %
%                                                                         %
%   Author  :   Gael de Oliveira                                          %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean Environment (a bit brutally, for sure!)
clear all; close all; clc; %#ok<CLALL>


% % Inputs are adimensionalized to references:
%       u_inf = 1     m/s       Freestream Speed
%       rho   = 1.225 kg/m3     Fluid Density (incomp.)
%       d     = 1     m         Diameter
Ct                  =  8/9;     % [-- ] Actuator Force Coefficient (C_F_a)
x_v                 = -0.2;     % [-- ] x-stance of lifting vortex pair (symmetry on x axis)
r_v                 =  0.4;     % [-- ] distance of lifting vortex pair to x axis (symmetry on x axis)
gamma_v             =  0.0;     % [-- ] strenght of each vortex in pair

% Post process
a          = 0.5 - 0.5 * sqrt(1 - Ct);
Se_over_Sa =      (1-a)/(1-2*a) ;
Re_over_Ra = sqrt((1-a)/(1-2*a));
% % Create Solver Object
DSD = kirikou_single_actuator_solver(Ct, x_v, r_v, gamma_v);

DSD.information_period  = 200;
DSD.theta_end           = 0.99*pi()/2;
DSD.max_RES_shape       = 1e-9;
DSD.max_RES_stretch     = 1e-9;
DSD.relax               = 0.010
DSD.dt                  = 0.010
DSD.n_stances           = 60;
% % Solve!
%DSD.preprocess_run_postprocess();

% % Prepocess Inputs
DSD.preprocess_inputs();
% % Discretize with Straight Wake Tube as Initial Guess
DSD.define_discretization_and_initial_geometry();
% % Define Reference Panel Strenghts
DSD.define_reference_panel_strenghts();
% % Now formulate Initial Guesses
DSD.formulate_initial_guesses();
% % Create Vortex Element Object (with state provided by initial guesses)
DSD.create_vortex_element_object();

%% Solution Process
% Now Start Solver
DSD.run_solver();
disp('Solver done!')
r_theo = DSD.y_start(1) * Re_over_Ra ;
r_calc = DSD.y_start(DSD.n_stances-1);
e_r    = (r_calc - r_theo) / r_theo  ;

disp(['r_theo = ' , num2str(r_theo )       ])
disp(['r_calc = ' , num2str(r_calc )       ])
disp(['e_r    = ' , num2str(e_r*100) , ' %'])



%% Postprocessing
% Now compute Machine Parameters
DSD.compute_machine_parameters();
% Generate Velocity Fields
DSD.generate_velocity_fields();

% % Enjoy results!
DSD.plot_velocity_field_with_streamlines();
print -dpng velocity_fields.png

% % Explore the fields of the DSD object, with a particular focus on the
% DSD.RES structure, which summarizes case results!
% Typical DSD.RES values (for example inputs):
%             r_a: 0.5000               % Actuator Radius
%           phi_a: -0.4444              % Actuator Force (Loading) Density
%             x_v: -0.2000              % x-stance of lifting vortices
%             r_v: 0.4000               % y-stance of lifting vortices
%         gamma_v: 0.4000               % circulation strenght of lifting vortices
%           C_F_a: 0.8889               % Actuator Force Coefficient
%           C_F_b: 0.2088               % Vortices Streamwise Force Coefficient (computed with generalized Kutta-Joukowski-Lagally theorem using numerical velocity)
%          Cp_num: -0.7348              % Numerical Cp, computed by integrating velocity field over actuator
%         Cp_theo: -0.7318              % Theoretical Cp, computed from C_F_b using the de Vries power coefficient law (as in my Torque 2016)
%            e_Cp: 0.0041               % Absolute Difference between Numerical Cp values
%          ua_num: 0.8266               % Average velocity on actuator plane (numerical)
%         ua_theo: 0.8233               % Average velocity on actuator plane (theorethical, given C_F_a and C_F_b)
%          a1_num: 0.1734               % Induction Factor on actuator plane (numerical)
%         a1_theo: 0.1767               % Induction Factor on actuator plane (theoretical)
%     RES_stretch: 9.9787e-13           % Solution Residual Filament Vorticity Evolution (RMS)
%       RES_shape: 8.8902e-13           % Solution Residual Wake Flow Tangency Condition (RMS)
%              VS: [1x1 constant_strenght_vortex_segment_2d]

%% Redefine Force Distribution
% Compute force on patches from force steps!
f1_actual_patches = - cumsum(DSD.delta_f_w .* sign(DSD.y_release));
% Reinterpolate through it with nearest neighbour interpolar
f1_actual = @(y) interp1(DSD.y_release, f1_actual_patches , y, 'next', 0);

%% Postprocessing (Compute Induction and Power Coefficient of Actuator) (New Method)

deltaR       = 0.000001;
DSD.x1       = 0;
%y_range1 = linspace(-(DSD.r1-eps), (DSD.r1-eps), 10000);
y_range1 = linspace(-(DSD.r1-deltaR), (DSD.r1-deltaR), 10000);
%y_range1 = linspace(-(DSD.r1-deltaR), (DSD.r1-deltaR), 1000000);
x_range1 = DSD.x1 * ones(size(y_range1));
% Compute wake effects on main extractor (all wakes!)
[~ , u_actuator1   , ~]   = induced_speed_on_many_points(DSD.VS, x_range1 , y_range1);
% Compute normal speed on actuator (only x component is needed,
% as we would take dot product of complete vector with actuator normal
% which is aligned with x axis!)
u_a1 = u_actuator1 + DSD.u_inf;
% Now compute induction factor
% a_actuator   = mean(u_inf-u_a_actuator(not(isnan(u_a_actuator))));
DSD.a1   = trapz(y_range1 , DSD.u_inf - u_a1) ./ (DSD.d1 *DSD.u_inf);
% Compute mass flow on surface
%m1   = trapz(y_range1, u_a1);
% And perform integral for power!
DSD.p1   = trapz(y_range1, f1_actual(y_range1) .* u_a1);
% Reach power coefficient!
DSD.Cp1  = DSD.p1 / (0.5*DSD.rho*DSD.S1*(DSD.u_inf^3));
% Display Results
Cp_num  = DSD.Cp1;
Cp_theo = -16/27;
e_Cp    = (Cp_num - Cp_theo) / Cp_theo;
disp(['Cp_num  = ' , num2str(Cp_num)])
disp(['Cp_theo = ' , num2str(Cp_theo)])
disp(['e_Cp    = ' , num2str(e_Cp)])


% Comparison of deltaR with distance between (elemental) vortex rings
% (smaller value is better!)
K = 128;
deltaR / (sqrt((DSD.x_start(2)-DSD.x_start(1)).^2 + (DSD.y_start(2)-DSD.y_start(1)).^2) / K)

%% Save stuff
save kirikou_example_case_1_axisymmetric_debug_alfa2.mat

%% Plot interesting stuff
figure(1985)
subplot(3,1,[1 2])
mu_range = [0, 0.2 , 0.4 , 0.6, 0.8, 0.95];
plot(DSD.x_mesh(1,:), interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh, DSD.x_mesh(1,:), 0           * (1/2)), '.-'); hold on
plot(DSD.x_mesh(1,:), interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh, DSD.x_mesh(1,:), mu_range(2) * (1/2)), '.-'); grid on
plot(DSD.x_mesh(1,:), interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh, DSD.x_mesh(1,:), mu_range(3) * (1/2)), '.-')
plot(DSD.x_mesh(1,:), interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh, DSD.x_mesh(1,:), mu_range(4) * (1/2)), '.-')
plot(DSD.x_mesh(1,:), interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh, DSD.x_mesh(1,:), mu_range(5) * (1/2)), '.-')
plot(DSD.x_mesh(1,:), interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh, DSD.x_mesh(1,:), mu_range(6) * (1/2)), '.-')

legend( 'Centerline', ...
        ['mu = ' , num2str(mu_range(2))], ...
        ['mu = ' , num2str(mu_range(3))], ...
        ['mu = ' , num2str(mu_range(4))], ...
        ['mu = ' , num2str(mu_range(5))], ...
        ['mu = ' , num2str(mu_range(6))])
    
axis([-1 , 1 , 1/3, 1]);
xlabel('x/D');
ylabel('u/U_{0}');
print -dpng axial_velocity_dependency_on_streamwise_stance.png

figure(1986)
mu_range = [0, 0.2 , 0.4 , 0.6, 0.8, 0.95];
plot(DSD.x_mesh(1,:), interp2(DSD.x_mesh, DSD.y_mesh, DSD.v_mesh, DSD.x_mesh(1,:), 0           * (1/2)), '.-'); hold on
plot(DSD.x_mesh(1,:), interp2(DSD.x_mesh, DSD.y_mesh, DSD.v_mesh, DSD.x_mesh(1,:), mu_range(2) * (1/2)), '.-'); grid on
plot(DSD.x_mesh(1,:), interp2(DSD.x_mesh, DSD.y_mesh, DSD.v_mesh, DSD.x_mesh(1,:), mu_range(3) * (1/2)), '.-')
plot(DSD.x_mesh(1,:), interp2(DSD.x_mesh, DSD.y_mesh, DSD.v_mesh, DSD.x_mesh(1,:), mu_range(4) * (1/2)), '.-')
plot(DSD.x_mesh(1,:), interp2(DSD.x_mesh, DSD.y_mesh, DSD.v_mesh, DSD.x_mesh(1,:), mu_range(5) * (1/2)), '.-')
plot(DSD.x_mesh(1,:), interp2(DSD.x_mesh, DSD.y_mesh, DSD.v_mesh, DSD.x_mesh(1,:), mu_range(6) * (1/2)), '.-')

legend( 'Centerline', ...
        ['mu = ' , num2str(mu_range(2))], ...
        ['mu = ' , num2str(mu_range(3))], ...
        ['mu = ' , num2str(mu_range(4))], ...
        ['mu = ' , num2str(mu_range(5))], ...
        ['mu = ' , num2str(mu_range(6))])
axis([-1, 1, 0, 0.3]);
xlabel('x/D');
ylabel('v/U_{0}');


figure(1987)
mu_range = [0, 0.2 , 0.4 , 0.6, 0.8, 0.95];
plot(DSD.x_mesh(1,:), atan2(interp2(DSD.x_mesh, DSD.y_mesh, DSD.v_mesh, DSD.x_mesh(1,:), 0           * (1/2)), interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh, DSD.x_mesh(1,:), 0           * (1/2))) * 180/pi(), '.-'); hold on
plot(DSD.x_mesh(1,:), atan2(interp2(DSD.x_mesh, DSD.y_mesh, DSD.v_mesh, DSD.x_mesh(1,:), mu_range(2) * (1/2)), interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh, DSD.x_mesh(1,:), mu_range(2) * (1/2))) * 180/pi(), '.-'); grid on
plot(DSD.x_mesh(1,:), atan2(interp2(DSD.x_mesh, DSD.y_mesh, DSD.v_mesh, DSD.x_mesh(1,:), mu_range(3) * (1/2)), interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh, DSD.x_mesh(1,:), mu_range(3) * (1/2))) * 180/pi(), '.-'); 
plot(DSD.x_mesh(1,:), atan2(interp2(DSD.x_mesh, DSD.y_mesh, DSD.v_mesh, DSD.x_mesh(1,:), mu_range(4) * (1/2)), interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh, DSD.x_mesh(1,:), mu_range(4) * (1/2))) * 180/pi(), '.-'); 
plot(DSD.x_mesh(1,:), atan2(interp2(DSD.x_mesh, DSD.y_mesh, DSD.v_mesh, DSD.x_mesh(1,:), mu_range(5) * (1/2)), interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh, DSD.x_mesh(1,:), mu_range(5) * (1/2))) * 180/pi(), '.-'); 
plot(DSD.x_mesh(1,:), atan2(interp2(DSD.x_mesh, DSD.y_mesh, DSD.v_mesh, DSD.x_mesh(1,:), mu_range(6) * (1/2)), interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh, DSD.x_mesh(1,:), mu_range(6) * (1/2))) * 180/pi(), '.-'); 


legend( 'Centerline', ...
        ['mu = ' , num2str(mu_range(2))], ...
        ['mu = ' , num2str(mu_range(3))], ...
        ['mu = ' , num2str(mu_range(4))], ...
        ['mu = ' , num2str(mu_range(5))], ...
        ['mu = ' , num2str(mu_range(6))])
    
xlabel('x/D');
ylabel('atan(v/u) [deg]');


figure(1988)
subplot(3,1,[1 2])
mu_range = [0, 0.2 , 0.4 , 0.6, 0.8, 0.95];
plot(DSD.x_mesh(1,:), 1 - interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh, DSD.x_mesh(1,:), 0           * (1/2)), '.-'); hold on
plot(DSD.x_mesh(1,:), 1 - interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh, DSD.x_mesh(1,:), mu_range(2) * (1/2)), '.-'); grid on
plot(DSD.x_mesh(1,:), 1 - interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh, DSD.x_mesh(1,:), mu_range(3) * (1/2)), '.-')
plot(DSD.x_mesh(1,:), 1 - interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh, DSD.x_mesh(1,:), mu_range(4) * (1/2)), '.-')
plot(DSD.x_mesh(1,:), 1 - interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh, DSD.x_mesh(1,:), mu_range(5) * (1/2)), '.-')
plot(DSD.x_mesh(1,:), 1 - interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh, DSD.x_mesh(1,:), mu_range(6) * (1/2)), '.-')

legend( 'Centerline', ...
        ['mu = ' , num2str(mu_range(2))], ...
        ['mu = ' , num2str(mu_range(3))], ...
        ['mu = ' , num2str(mu_range(4))], ...
        ['mu = ' , num2str(mu_range(5))], ...
        ['mu = ' , num2str(mu_range(6))])
    
xlabel('x/D');
ylabel('a = 1 - u/U_{0}');
axis([-2 2 0 2/3])
print -dpng a_dependency_on_streamwise_stance.png

figure(1989)
subplot(3,1,[1 2])
mu_range = [0, 0.2 , 0.4 , 0.6, 0.8, 0.95];
plot(DSD.y_mesh(:,1), 1 - interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh,-0.08 , DSD.y_mesh(:,1)), '.-'); hold on
plot(DSD.y_mesh(:,1), 1 - interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh,-0.04 , DSD.y_mesh(:,1)), '.-'); grid on
plot(DSD.y_mesh(:,1), 1 - interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh, 0    , DSD.y_mesh(:,1)), '.-');
plot(DSD.y_mesh(:,1), 1 - interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh, 0.04 , DSD.y_mesh(:,1)), '.-');
plot(DSD.y_mesh(:,1), 1 - interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh, 0.08 , DSD.y_mesh(:,1)), '.-');

legend('x/D=-0.08', 'x/D=-0.04', 'x/D=-0', 'x/D= 0.04', 'x/D=-0.08')
    
xlabel('r/D');
ylabel('a = 1 - u/U_{0}');

figure(2018)
subplot(3,1,[1 2])
plot(DSD.y_mesh(:,1), 1 - interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh, 0    , DSD.y_mesh(:,1)), '.-'); hold on;
plot([-1 1], [1 1]/3, 'k--')
axis([-0.5, 0.5, 0.2, 0.5]); grid on
xlabel('r/D');
ylabel('u/U_{0}');
legend('Numerical solution', 'Average=BEM')
print -dpng axial_velocity_along_radius.png

%% And plot more interesting stuff
D = 1/2;
coning_angle_range = [0 , 5, 10] * pi/180;
y_range            = linspace(0,1/2,100);

x_range1 = -y_range * tan(coning_angle_range(1)) + tan(coning_angle_range(1))*D;
x_range2 = -y_range * tan(coning_angle_range(2)) + tan(coning_angle_range(2))*D;
x_range3 = -y_range * tan(coning_angle_range(3)) + tan(coning_angle_range(3))*D;

figure(2019)
subplot(3,1,[1 2])
n_pps = length(DSD.x_start)/2;
h1 = plot(x_range1 , y_range); hold on; 
h2 = plot(x_range2 , y_range); 
h3 = plot(x_range3 , y_range); 
plot(DSD.x_start(1:n_pps), DSD.y_start(1:n_pps), 'ko-'); hold on;
plot(DSD.x_start(n_pps+(1:n_pps)), DSD.y_start(n_pps+(1:n_pps)), 'ko-');
plot(x_range1 , - y_range, 'Color', h1.Color); 
plot(x_range2 , - y_range, 'Color', h2.Color); 
plot(x_range3 , - y_range, 'Color', h3.Color); 

grid on;
axis([-1 2 -0.8 0.8])
legend( 'coning angle = 0', ...
       ['coning angle = '  , num2str(coning_angle_range(2)*180/pi), 'deg'] , ...
       ['coning angle = '  , num2str(coning_angle_range(3)*180/pi), 'deg'] , ...
        'wake edges');
print -dpng effect_of_coning_on_induction_sampling_points.png

figure(1990)
subplot(3,1,[1 2])
plot(y_range, 1 - interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh,x_range1, y_range), '.-'); hold on
plot(y_range, 1 - interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh,x_range2, y_range), '.-'); grid on
plot(y_range, 1 - interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh,x_range3, y_range), '.-');

legend( 'coning angle = -0', ...
       ['coning angle = '  , num2str(coning_angle_range(2)*180/pi)] , ...
       ['coning angle = '  , num2str(coning_angle_range(3)*180/pi)] );
axis([0 , 0.5 , 0.25, 0.5])
xlabel('y/d');
ylabel('a = 1 - u/U_{0}');
print -dpng effect_of_coning_on_induction.png

U0      = 1;
lambda  = 8;
% phi   = lambda ./ sqrt((1-a).^2 + lambda.^2);

a1 = 1 - interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh,x_range1, y_range);
a2 = 1 - interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh,x_range2, y_range);
a3 = 1 - interp2(DSD.x_mesh, DSD.y_mesh, DSD.u_mesh,x_range3, y_range);


phi1 = atan((1-a1) ./ (lambda*y_range).^2 );
phi2 = atan((1-a2) ./ (lambda*y_range).^2 );
phi3 = atan((1-a3) ./ (lambda*y_range).^2 );


figure(1991)
subplot(3,1,[1 2])
plot(y_range, phi1 * 180/pi(), '.-'); hold on
plot(y_range, phi2 * 180/pi(), '.-'); grid on
plot(y_range, phi3 * 180/pi(), '.-');
legend( 'coning angle = -0', ...
       ['coning angle = '  , num2str(coning_angle_range(2)*180/pi())] , ...
       ['coning angle = '  , num2str(coning_angle_range(3)*180/pi())] );
axis([0 0.5 -5 30])
xlabel('y/d');
ylabel('\phi = atan((1-a)/lambda) [deg]');
title(['Effect of conning on induction at rotor (\lambda = ' , num2str(lambda), ') : inflow angle [deg]'])
print -dpng effect_of_coning_on_induction_angle.png

figure(1992)
subplot(3,1,[1 2])
plot(y_range, (phi1 - phi1) * 180/pi(), '.-'); hold on
plot(y_range, (phi2 - phi1) * 180/pi(), '.-'); grid on
plot(y_range, (phi3 - phi1) * 180/pi(), '.-');
legend( 'coning angle = -0', ...
       ['coning angle = '  , num2str(coning_angle_range(2)*180/pi()) , ' (positive = tip upwind)'] , ...
       ['coning angle = '  , num2str(coning_angle_range(3)*180/pi()) , ' (positive = tip upwind)'] );
axis([0 0.5 -2.0 0.5])
xlabel('y/d');
ylabel('\phi - \phi_{no coning}');
title(['Effect of conning on induction at rotor (\lambda = ' , num2str(lambda), ') : inflow angle offset [deg]'])
print -dpng effect_of_coning_on_induction_angle_offset.png

figure(1993)
subplot(3,1,[1 2])
plot(y_range, 2 * pi() * sin(phi1 - phi1), '.-'); hold on
plot(y_range, 2 * pi() * sin(phi2 - phi1), '.-'); grid on
plot(y_range, 2 * pi() * sin(phi3 - phi1), '.-');
legend( 'coning angle = -0', ...
       ['coning angle = '  , num2str(coning_angle_range(2)*180/pi()) , ' (positive = tip upwind)'] , ...
       ['coning angle = '  , num2str(coning_angle_range(3)*180/pi()) , ' (positive = tip upwind)'] );
axis([0 0.5 -0.25 0.05])
xlabel('y/d');
ylabel('C_l offset (2\pi sin \alpha');
title(['Effect of conning on induction at rotor (\lambda = ' , num2str(lambda), ') : lift coefficient offset'])
print -dpng effect_of_coning_on_induction_angle_cl_offset.png





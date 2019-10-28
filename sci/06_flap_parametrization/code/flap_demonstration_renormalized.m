


% Close plots
close all

% Generate Coordinates
[tx, tz] = SD.generate_coordinates(400, x);

% Set Hinge chordwise position
SDD.x_hinge = 0.6;
% Set knee parameter
SDD.flap_knee = 0.1;
% Set flap angle
flap_angle = 20;

% Apply transformation Down
SDD.flap_angle = flap_angle;       % Set flap angle
[tsd, tnd] = SDD.T_curvilinear_transformation(tx, tz);

% Apply transformation Up
SDD.flap_angle = -flap_angle;       % Set flap angle
[tsu, tnu] = SDD.T_curvilinear_transformation(tx, tz);

% And plot the two airfoils!
plot(tx, tz, tsd, tnd, tsu, tnu); 
grid on; axis([0 , 1 , -0.4, 0.4])
hold on

% Add a circle to highlight T.E position
r = 1-SDD.x_hinge;                      % Circle Radius
t = 0.25 * linspace(-pi(), pi());       % Angle Range

plot(r*cos(t) + (1-r), r*sin(t));       % Plot Circle 
plot( 1-r , 0 , 'o');                   % Mark Hinge position

% Beautify Plot
legend('Original Airfoil' , 'Flapped Down \delta=20°', 'Flapped Up \delta=-20°', 'T.E. Locii', 'Hinge Location' , 'Location', 'South')
title('Renormalized quasi-affine transformation - Merci Lamé (1859)')
xlabel('s')
ylabel('\eta')

% Now store plot to file
set(gcf, 'PaperType', 'A5')
orient landscape
print('-dpdf', 'renormalized_airfoil_flapping.pdf')
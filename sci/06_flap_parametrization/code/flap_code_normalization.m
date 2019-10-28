%clear all; 
close all; clc
% Inputs: hinge location and smoothing parameter
x_hinge     = 0.5;
z_hinge     = 0;
flap_knee   = 0.001;
flap_angle  = 15; 

% Create Shape Dynamizer object!
SDD = shape_dynamizer();
SDD.x_hinge     = x_hinge;
SDD.z_hinge     = z_hinge;
SDD.flap_knee   = flap_knee;
SDD.flap_angle  = flap_angle;

% Make a reference hingeline vector
tx = linspace(0,1,1000);
tz = ones(size(tx)) * (0);

% Now transform the hingeline
[ts,tn] = SDD.T_curvilinear_transformation(tx, tz);
tz2 = SDD.hingeline_deformation(tx);

% And plot
plot(tx, tz, tx, tz2, '--', ts, tn);
grid on
legend('Original' , 'Deformed', 'Deformed and Renormalized' )
axis([-0.1 1.1 -0.6 0.6])


disp('Original Lenght and Deformed length')
sum(sqrt(diff(tx).^2 + (diff(tz)).^2))
sum(sqrt(diff(ts).^2 + (diff(tn)).^2))
disp('Effective hingeline angle:')
atan2(tn(end) - tz(1), ts(end)-x_hinge) * 180 / pi()
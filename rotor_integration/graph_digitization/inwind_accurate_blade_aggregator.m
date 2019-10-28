% Load data grabed by hand from grabbit
load inwind_chord.mat
load inwind_thickness.mat
load inwind_twist_deg.mat

inwind_RWT = struct();

inwind_RWT.lambda_design = 7.5;
inwind_RWT.R             = 89.166;
inwind_RWT.R_hub         = 5.6/2;

inwind_RWT.r_chord       = [0; inwind_chord(2:(end),1); inwind_RWT.R] + inwind_RWT.R_hub;
inwind_RWT.val_chord     = [5.3805; inwind_chord(2:(end),2); 0];

inwind_RWT.r_thickness   = [0 ; inwind_thickness(2:end,1); inwind_RWT.R] + inwind_RWT.R_hub;
inwind_RWT.val_thickness = [100; 100; inwind_thickness(3:end,2); inwind_thickness(end, 2)];

inwind_RWT.r_twist_deg   = [0 ; inwind_twist_deg(1:end, 1); inwind_RWT.R] + inwind_RWT.R_hub;
inwind_RWT.val_twist_deg = [inwind_twist_deg(1, 2) ; inwind_twist_deg(1:end, 2); 0];
inwind_RWT.R             = 89.166 + inwind_RWT.R_hub;

save('inwind_RWT.mat', 'inwind_RWT');

figure(1)
subplot(311)
plot(inwind_RWT.r_chord, inwind_RWT.val_chord)
xlabel('Radius (m)'); ylabel('Chord (m)'); grid on
subplot(312)
plot(inwind_RWT.r_thickness, inwind_RWT.val_thickness)
xlabel('Radius (m)'); ylabel('Thickness (%)'); grid on
subplot(313)
plot(inwind_RWT.r_twist_deg, inwind_RWT.val_twist_deg)
xlabel('Radius (m)'); ylabel('Twist (deg)'); grid on

figure(2)
subplot(311)
plot(inwind_RWT.r_chord / inwind_RWT.R, inwind_RWT.val_chord / inwind_RWT.R * 100)
xlabel('Radius (r/R)'); ylabel('c/R (%)'); grid on
subplot(312)
plot(inwind_RWT.r_thickness / inwind_RWT.R, inwind_RWT.val_thickness)
xlabel('Radius (r/R)'); ylabel('t/c (%)'); grid on
subplot(313)
plot(inwind_RWT.r_twist_deg / inwind_RWT.R, inwind_RWT.val_twist_deg)
xlabel('Radius (r/R)'); ylabel('Twist (deg)'); grid on







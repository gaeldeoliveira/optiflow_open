% Make 

% Data from Rahim' paper
% stdev of AOA, at different radial positions [30,40,50,60,70,80] m in each
% column
% and differenet cases in each line ( 3, 5, 9, 30, 33)
data_Rahim=[0.7 0.5 0.45 0.4 0.4 0.4;0.6 0.6 0.65 0.7 0.75 0.75; 3.25 2.65 2.3 1.75 1.65 1.6; 3.3 2.7 2.2 1.95 1.85 1.8;3.9 3.3 2.7 2.35 2.25 2];
% Now extract components, as per Torque paper's graph's (and replicate at
% edges (0 and R=91.966m) to avoid numerical interpolation issues)
r_range             = [0-eps(1) 30 40 50 60 70 80 (89.166+5.6/2)+eps(100)];
std_alpha_deg_yaw   = [data_Rahim(1,1), data_Rahim(1,:), data_Rahim(1,end)];
std_alpha_deg_shear = [data_Rahim(2,1), data_Rahim(2,:), data_Rahim(2,end)];
std_alpha_deg_turb  = [data_Rahim(3,1), data_Rahim(3,:), data_Rahim(3,end)];

% Combine perturbations as if all were normal distributions
std_alpha_deg_total = sqrt(std_alpha_deg_yaw.^2 + std_alpha_deg_shear.^2 + std_alpha_deg_turb.^2);

% Bundle all this into a structure
std_distribution = struct();
std_distribution.r_range             = r_range           ;
std_distribution.mu_range            = r_range/(89.166+5.6/2);
std_distribution.std_alpha_deg_total = std_alpha_deg_total;
std_distribution.std_alpha_deg_yaw   = std_alpha_deg_yaw  ;
std_distribution.std_alpha_deg_shear = std_alpha_deg_shear;
std_distribution.std_alpha_deg_turb  = std_alpha_deg_turb ;
std_distribution.source              = 'RE paper of Rahim and RIcardo, as reused in our Torque 2018 paper';  %#ok<STRNU>
% And save
save('std_distribution.mat', 'std_distribution');

std_alpha_deg_turb_scaled  = std_alpha_deg_turb * 8.4 / 7.5;
std_alpha_deg_total_scaled = sqrt(std_alpha_deg_yaw.^2 + std_alpha_deg_shear.^2 + std_alpha_deg_turb_scaled.^2);

% Bundle all this into a structure
std_distribution = struct();
std_distribution.r_range             = r_range           ;
std_distribution.mu_range            = r_range/(89.166+5.6/2);
std_distribution.std_alpha_deg_total = std_alpha_deg_total_scaled;
std_distribution.std_alpha_deg_yaw   = std_alpha_deg_yaw  ;
std_distribution.std_alpha_deg_shear = std_alpha_deg_shear;
std_distribution.std_alpha_deg_turb  = std_alpha_deg_turb_scaled ;
std_distribution.source              = 'RE paper of Rahim and RIcardo, as reused in our Torque 2018 paper, TI upscaled with TSR from 8.4 (Rahim) to 7.5 (our case).';
% And save
save('std_distribution_scaled.mat', 'std_distribution');
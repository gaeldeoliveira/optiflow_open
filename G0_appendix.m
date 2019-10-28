% Draft for code to be added at end of runfiles
m_str = mfilename();
save([m_str, '-results-', datestr(now, 'yyyymmddTHHMMSS')], 'results');
save([m_str, '-fullcase-', datestr(now, 'yyyymmddTHHMMSS')]);


% Draft for code to generate angle of attach range
% Wind and Kite Condition Inputs
u_wind    =  10;                        % [m/s ] wind velocity
t_u       = 0.1;                        % [m/s ] 10 percent relative wind turbulence (RMS)
u_kite    =  30;                        % [m/s ] kite velocity
c_kite    =  0.3;                       % [m   ] representative chord of kite wing
nu_air    =  1.48e-5;                   % [m2/s] kinematic viscosity of air (at 15C, take 1.57 for 25C)
cl_design =  1.2;

sigma_u = u_wind * t_u;                 % [m/s] RMS of wind

phi_avg = atan2(u_wind, u_kite);        % [rad] average inflow angle
% Very crude estimate of RMS of inflow angle, we have a much better way to
% do this, which considers the fact that this distribution is not gaussian
phi_RMS  = 0.5*(atan2(u_wind+sigma_u, u_kite) - atan2(u_wind-sigma_u, u_kite));

% And just as crude estimate of RMS of lift
cl_RMS   = 2*pi()*phi_RMS;
eq_stall_margin = (cl_RMS*3/2) / (2*pi()) * (180/pi); % (at 3/2 sigma)

% Define normal distribution functions
norm_pdf = @(x) 1/sqrt(2*pi()) * exp(-(x/2).^2);
norm_cdf = @(x) 0.5*(1+erf(x/sqrt(2)));

% Define vector of standard deviation steps
sigma_vector     = (1/4) * [-6 -4 -2 0 2 4 6] ;
% Compute Weight Vector (by hand, for sigma_vector = [-6 -4 -2 0 2 4 6]/4)
cl_weight_vector = [norm_cdf(-5/4)-norm_cdf(-5  ) , ...
                    norm_cdf(-3/4)-norm_cdf(-5/4) , ...
                    norm_cdf(-1/4)-norm_cdf(-3/4) , ...
                    norm_cdf( 1/4)-norm_cdf(-1/4) , ...
                    norm_cdf( 3/4)-norm_cdf( 1/4) , ...
                    norm_cdf( 5/4)-norm_cdf( 3/4) , ...
                    norm_cdf( 5  )-norm_cdf( 5/4)];
% Define vectors of lift coefficient steps
cl_sigma_vector  = cl_RMS * sigma_vector;               % centered on 0
cl_point_vector   = cl_RMS * sigma_vector + cl_design;   % centered on cl_design


disp('-----------------------------------------------------------------------')
disp(['Kite Speed      : ' , num2str(u_kite)           , ' m/s'])
disp(['Wind Speed      : ' , num2str(u_wind)           , ' m/s'])
disp(['Turb. Intensity : ' , num2str(t_u*100)          , ' percent'])
disp('-----------------------------------------------------------------------')
disp(['Inflow Angle Avg: ' , num2str(phi_avg*180/pi()) , ' deg'])
disp(['Inflow Angle RMS: ' , num2str(phi_RMS*180/pi()) , ' deg'])
disp(['Linear lift  RMS: ' , num2str(cl_RMS)           , ' adim.'])
disp('-----------------------------------------------------------------------')
disp(['Eq. stall margin: ' , num2str(eq_stall_margin ) , ' deg (3/2 sigma)'])
disp('Evaluation Ranges, [Weights ; StDev (RMS) Steps ; Cl Steps; Cl Points] : ')
disp([weight_vector ; phi_sigma_vector ; cl_sigma_vector; cl_point_vector])
disp('-----------------------------------------------------------------------')


                
                





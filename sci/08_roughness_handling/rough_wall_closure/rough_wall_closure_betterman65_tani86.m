function [H12, Re_theta_S, Re_theta_R] = rough_wall_closure_betterman65_tani86(cf_input, Pi_wake, Re_h, lambda)

% General constants
kappa       =  0.41   ;    % [adim.] Von Karmann constant
C           =  5.0    ;    % [adim.] Smooth wall intersect

% Map skin friction into 
z           = kappa * sqrt(2/cf_input);

% % Compute intermediate variables
alfa        = 59/60 + Pi_wake;
beta        = 8437/4200 + 667/210 * Pi_wake + 52/35*Pi_wake.^2;
s           = z - 2*Pi_wake - kappa * C;

% % Compute K constant of roughness shift function 
% lambda = linspace(1, 20);
% Compute function branches separately
f_low_lambda  =  17.35 * (1.634 * log10(lambda) - 1);
f_high_lambda = - 5.95 * (1.103 * log10(lambda) - 1);
% Define branch threshold (get it like this, lambda_threshold = fzero(@(lambda_threshold) f_low_lambda(lambda_threshold) - f_high_lambda(lambda_threshold), 4.6491972))
lambda_threshold = 4.6491972;
% Allocate solutions
f_lambda      = zeros(size(lambda));
% Join the branches
f_lambda(lambda>lambda_threshold)  = f_high_lambda(lambda>lambda_threshold);
f_lambda(lambda<=lambda_threshold) = f_low_lambda(lambda<=lambda_threshold);

% % Compute roughness shift function
DeltaU_plus = (1/kappa) .* log((kappa./z).*Re_h) + f_lambda;
% Fudge points where 
DeltaU_plus(DeltaU_plus<0) = 0;

% % Test plot (for roughness
% Set lambda as:
%lambda = linspace(1, 20);
% (re)Run necessary lines and plot:
%plot(lambda, f_lambda); hold on; plot(lambda, DeltaU_plus ); grid on; legend('f of \lambda', '\Delta U / u_\tau');

% % Compute Shape Factor
H12         = alfa ./ (alfa - beta ./ z);
% % Compute Re_theta for smooth case (with this z, PIwake combination) 
Re_theta_S  = (alfa - beta ./ z) .* (1/kappa) .* exp(s                      );
% % Compute Re_theta for rough case (with this z, PIwake combination) 
Re_theta_R  = (alfa - beta ./ z) .* (1/kappa) .* exp(s + kappa * DeltaU_plus);

end
function [H12_S, Re_theta_S, cf_S, exitflag_S, H12_R, Re_theta_R, cf_R, exitflag_R] = invert_rough_wall_closure_betterman65_tani86(H12_target, Re_theta_target, cf0, Pi_wake0, Re_h, lambda, lsq_options)

% % Multiplex starting point
x0(1) = cf0;
x0(2) = Pi_wake0;
% % Multiplex target point
y_target(1)     = H12_target            ;
y_target(2)     = Re_theta_target / 1000;

% Make scalar residual function
% res_S = @(x) sqrt(sum(( y_target - rough_wall_closure_betterman65_tani86_wrapper_S(x, Re_h, lambda) ).^2));
% res_R = @(x) sqrt(sum(( y_target - rough_wall_closure_betterman65_tani86_wrapper_R(x, Re_h, lambda) ).^2));
% Solve as minimization problem
% [x_smooth,fval_smooth,exitflag_smooth,output_smooth,lambda_smooth,grad_smooth,hessian_smooth] = fmincon(res_S ,x0      , [],[],[],[],[0 0],[1, 20]);
% [x_rough ,fval_rough ,exitflag_rough ,output_rough ,lambda_rough ,grad_rough ,hessian_rough ] = fmincon(res_R ,x_smooth, [],[],[],[],[0 0],[1, 20]);

% Make pointwise (vector) residual function
res_S_LM = @(x) real(y_target - rough_wall_closure_betterman65_tani86_wrapper_S(real(x), Re_h, lambda));
res_R_LM = @(x) real(y_target - rough_wall_closure_betterman65_tani86_wrapper_R(real(x), Re_h, lambda));
% Solve as curve-fitting problem (Levenberg-Marquardt problem)
if isempty(lsq_options)
    [x_smooth,~,~,exitflag_S,~] = lsqnonlin(res_S_LM,x0);
    [x_rough ,~,~,exitflag_R,~] = lsqnonlin(res_R_LM,x0);
else
    [x_smooth,~,~,exitflag_S,~] = lsqnonlin(res_S_LM,x0, [], [], lsq_options);
    [x_rough ,~,~,exitflag_R,~] = lsqnonlin(res_R_LM,x0, [], [], lsq_options);
end


[H12_S, Re_theta_S, ~] = rough_wall_closure_betterman65_tani86(x_smooth(1), x_smooth(2), Re_h, lambda);
[H12_R, ~, Re_theta_R] = rough_wall_closure_betterman65_tani86(x_rough( 1), x_rough( 2), Re_h, lambda);

cf_S = x_smooth(1);
cf_R = x_rough( 1);

end
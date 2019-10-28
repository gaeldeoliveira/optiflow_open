%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%       Add k-type roughness to TBL treatment of Rfoil                    %
%       2nd approach, all fully based on Itiro Tani (boss of the bosses!) %
%           Gael de Oliveira                                              %
%           May 2018                                                      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs definning similarity conditions of specific point of interest
cf_input    = 0.001;       % [adim.] Skin friction coefficient of current boundary layer
Pi_wake     = 5    ;       % [adim.] Factor multiplying universal wake function (based on Coles rationale, but we use the lewkowicz 1982 polynomial instead, after Itiro Tani)
Re_h        = 50   ;       % [adim.] Reynolds number of roughness height (based no u_edge, not u_tau. We have Re_h = U_edge * h / nu , where h is height of roughess elements)
lambda      = 10   ;       % [adim.] Ratio (D/h) between roughness pitch (D) and height (d) (pas/hauteur, as defined in Betterman 1965, IJHMT vol. 9 pp. 153-165)

% % Grant access to rough wall closure helper functions
addpath rough_wall_closure/

% % Compute closure from velocity profile
[H12, Re_theta_S, Re_theta_R] = rough_wall_closure_betterman65_tani86(cf_input, Pi_wake, Re_h, lambda);

% % Display Results
disp(['Rough  >> H12 = ' , num2str(H12, 3), '     Re_theta = ', num2str(Re_theta_R, 5)]);
disp(['Smooth >> H12 = ' , num2str(H12, 3), '     Re_theta = ', num2str(Re_theta_S, 5)]);


% % Compare with current closure!
addpath dev/0_closure_relations/
% % Compare with current closure!
[ cf_reference ] = cft_seq( H12, Re_theta_S, 0);
[ cf_smooth_at_rough  ] = cft_seq( H12, Re_theta_R, 0);
% % Compute differences!
delta_cf_smooth    = cf_reference - cf_input;
delta_cf_roughness = cf_reference - cf_smooth_at_rough ;


% % Now move on to the computation of Cf's at prescribed H12 and Re_theta
% pairs (in other words, find 

% Define target point
H12_target      =    1.5;
Re_theta_target = 8000  ;
% % Set starting point
cf0 = cft_seq( H12_target, Re_theta_target, 0);
Pi_wake0 = 1;
% % Invert closure
% [H12_S, Re_theta_S, cf_S, H12_R, Re_theta_R, cf_R] = invert_rough_wall_closure_betterman65_tani86(H12_target, Re_theta_target, cf0, Pi_wake0, Re_h, lambda, []);
[H12_S, Re_theta_S, cf_S, exitflag_S, H12_R, Re_theta_R, cf_R, exitflag_R] = invert_rough_wall_closure_betterman65_tani86(H12_target, Re_theta_target, cf0, Pi_wake0, Re_h, lambda, []);
%lsq_options = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt', 'MaxIterations', 1000)
%[H12_S, Re_theta_S, cf_S, H12_R, Re_theta_R, cf_R] = invert_rough_wall_closure_betterman65_tani86(H12_target, Re_theta_target, cf0, Pi_wake0, Re_h, lambda, lsq_options );
% % Display results
disp(['H12_S = ', num2str(H12_S), '    Re_theta_S = ', num2str(Re_theta_S), '    Cf_S = ' , num2str(cf_S), '    flag_S = ' , num2str(exitflag_S)]);
disp(['H12_R = ', num2str(H12_R), '    Re_theta_R = ', num2str(Re_theta_R), '    Cf_R = ' , num2str(cf_R), '    flag_R = ' , num2str(exitflag_R)]);

(cf_R - cf_S) / cf_S


% % % % Now make a mesh of target points to generate numerical closure 
% % Define closure generation ranges
lambda_target         = 10;
Re_h_target_range     = [0 , 50, 100];
Re_theta_target_range = [800, 1600, 3200, 6400, 12800];
H12_target_range      = [1.4, 1.6, 1.8, 2.0, 2.2, 2.4];
% % And associated initial guesses for solver (subordinated to H12)
Pi_wake0_target_range = [1  , 1  , 2  , 2  , 6.5, 9.9];

% % Allocate nd-meshes (grids) for results (including consistency check data)
mesh_H12_S      = zeros(length(Re_theta_target_range), length(H12_target_range), length(Re_h_target_range));
mesh_H12_R      = zeros(length(Re_theta_target_range), length(H12_target_range), length(Re_h_target_range));
mesh_Re_theta_S = zeros(length(Re_theta_target_range), length(H12_target_range), length(Re_h_target_range));
mesh_Re_theta_R = zeros(length(Re_theta_target_range), length(H12_target_range), length(Re_h_target_range));
mesh_Cf_S       = zeros(length(Re_theta_target_range), length(H12_target_range), length(Re_h_target_range));
mesh_Cf_R       = zeros(length(Re_theta_target_range), length(H12_target_range), length(Re_h_target_range));
mesh_flag_S     = zeros(length(Re_theta_target_range), length(H12_target_range), length(Re_h_target_range));
mesh_flag_R     = zeros(length(Re_theta_target_range), length(H12_target_range), length(Re_h_target_range));
mesh_eps_Cf     = zeros(length(Re_theta_target_range), length(H12_target_range), length(Re_h_target_range));
mesh_Re_h       = zeros(length(Re_theta_target_range), length(H12_target_range), length(Re_h_target_range));


for n_H12_target = 1:length(H12_target_range)
    for n_Re_theta_target = 1:length(Re_theta_target_range)
        for n_Re_h_target = 1:length(Re_h_target_range)
            % % Get target conditions of current point
            Re_theta_target = Re_theta_target_range(n_Re_theta_target);
            H12_target      = H12_target_range(     n_H12_target);
            Re_h_target     = Re_h_target_range(    n_Re_h_target);
            
            % % Set initial guess for current target point
            cf0             = cft_seq( H12_target, Re_theta_target, 0);
            Pi_wake0        = Pi_wake0_target_range(n_H12_target);
            
            % % Invert closure for current target point
            [H12_S, Re_theta_S, cf_S, exitflag_S, ...
                H12_R, Re_theta_R, cf_R, exitflag_R] = ...
                invert_rough_wall_closure_betterman65_tani86( ...
                H12_target, Re_theta_target, cf0, Pi_wake0, Re_h_target, lambda_target, []);
            
            % % Compute effect on Skin friction ratio
            eps_Cf = (cf_R - cf_S) / cf_S;
            
            % % Display results
            disp(['H12_S = ', num2str(H12_S), '    Re_theta_S = ', num2str(Re_theta_S), '    Cf_S = ' , num2str(cf_S), '    flag_S = ' , num2str(exitflag_S)]);
            disp(['H12_R = ', num2str(H12_R), '    Re_theta_R = ', num2str(Re_theta_R), '    Cf_R = ' , num2str(cf_R), '    flag_R = ' , num2str(exitflag_R)]);
            
            % % Fill current point results into mesh
            mesh_H12_S(     n_Re_theta_target, n_H12_target, n_Re_h_target) = H12_S      ;
            mesh_H12_R(     n_Re_theta_target, n_H12_target, n_Re_h_target) = H12_R      ;
            mesh_Re_theta_S(n_Re_theta_target, n_H12_target, n_Re_h_target) = Re_theta_S ;
            mesh_Re_theta_R(n_Re_theta_target, n_H12_target, n_Re_h_target) = Re_theta_R ;
            mesh_Cf_S(      n_Re_theta_target, n_H12_target, n_Re_h_target) = cf_S       ;
            mesh_Cf_R(      n_Re_theta_target, n_H12_target, n_Re_h_target) = cf_R       ;
            mesh_flag_S(    n_Re_theta_target, n_H12_target, n_Re_h_target) = exitflag_S ;
            mesh_flag_R(    n_Re_theta_target, n_H12_target, n_Re_h_target) = exitflag_R ;
            mesh_eps_Cf(    n_Re_theta_target, n_H12_target, n_Re_h_target) = eps_Cf     ;
            mesh_Re_h(      n_Re_theta_target, n_H12_target, n_Re_h_target) = Re_h_target;
        end
    end
end

% Checked:
% For Re_h_target = 50 and lambda_target = 10 :
%   n = 1   H12 = 1.4   OK (for Re_theta = 800, 1600, 3200, 6400, 12800)
%   n = 2   H12 = 1.6   OK (for Re_theta = 800, 1600, 3200, 6400, 12800)
%   n = 3   H12 = 1.8   OK (for Re_theta = 800, 1600, 3200, 6400, 12800)
%   n = 4   H12 = 2.0   OK (for Re_theta = 800, 1600, 3200, 6400, 12800)
%   n = 5   H12 = 2.2   OK (for Re_theta = 800, 1600, 3200, 6400, 12800) (after calibration of starting condition, Pi_wake0 = 6.5)
%   n = 6   H12 = 2.2   OK (for Re_theta = 800, 1600, 3200, 6400, 12800) (after calibration of starting condition, Pi_wake0 = 9.9)
%
% For Re_h_target =100 and lambda_target = 10 :
%   n = 1   H12 = 1.4   OK (for Re_theta = 800, 1600, 3200, 6400, 12800)
%   n = 2   H12 = 1.6   OK (for Re_theta = 800, 1600, 3200, 6400, 12800)
%   n = 3   H12 = 1.8   OK (for Re_theta = 800, 1600, 3200, 6400, 12800)
%   n = 4   H12 = 2.0   OK (for Re_theta = 800, 1600, 3200, 6400, 12800)
%   n = 5   H12 = 2.2   OK (for Re_theta = 800, 1600, 3200, 6400, 12800) (after calibration of starting condition, Pi_wake0 = 6.5)
%   n = 6   H12 = 2.2   OK (for Re_theta = 800, 1600, 3200, 6400, 12800) (after calibration of starting condition, Pi_wake0 = 9.9)
%
% For Re_h_target =  0 and lambda_target = 10 :
%   n = 1   H12 = 1.4   OK (for Re_theta = 800, 1600, 3200, 6400, 12800)
%   n = 2   H12 = 1.6   OK (for Re_theta = 800, 1600, 3200, 6400, 12800)
%   n = 3   H12 = 1.8   OK (for Re_theta = 800, 1600, 3200, 6400, 12800)
%   n = 4   H12 = 2.0   OK (for Re_theta = 800, 1600, 3200, 6400, 12800)
%   n = 5   H12 = 2.2   OK (for Re_theta = 800, 1600, 3200, 6400, 12800)
%   n = 6   H12 = 2.2   OK (for Re_theta = 800, 1600, 3200, 6400, 12800)



% Johnny Haliday - Je Lis

% % Retrofit variable names to reuse Fortran exportation script that was
% written for Cei closure used in plasmas.
% Also, replicate solutions for hk = 1.4 to make them usable until hk = 1

% Match independent variable ranges
hk_range        = H12_target_range;
rt_range        = Re_theta_target_range;
treh_range      = Re_h_target_range;

% Match independent variable grids
hk_grid         = mesh_H12_S;
rt_grid         = mesh_Re_theta_S;
treh_grid       = mesh_Re_h;
% Match   dependent variable grid
cec_grid        = mesh_eps_Cf;

% % Remake grids to avoid loss of consistency due to residual solver errors
[hk_grid_b , rt_grid_b , t_theta_p_grid_b] = meshgrid(hk_range, rt_range, treh_range);
disp(['Global HK inversion error norm: ',  num2str(sqrt(sum(sum(sum((hk_grid - hk_grid_b).^2)))))]);
disp(['Global RT inversion error norm: ',  num2str(sqrt(sum(sum(sum((rt_grid - rt_grid_b).^2)))))]);
disp(['Replace (HK,RT )grids with remeshed ones']); %#ok<NBRAK>
hk_grid = hk_grid_b;
rt_grid = rt_grid_b;

% % Export to fortran code (float->symbolic->fortran string->code file)

% % Recast grid into table for human readable fortran code
cec_table = zeros(size(cec_grid, 2), size(cec_grid,1), size(cec_grid, 3));
for i=1:size(cec_table,1)
    for j=1:size(cec_table,2)
        for k=1:size(cec_table,3)
            cec_table(i,j,k) = real(cec_grid(j,i,k));
            if abs(imag(cec_grid(j,i,k))) > 0
                disp('WARNING: Imaginary Components ?(-1)^-3? truncated!')
            end
        end
    end
end

% % Start by making a carriage return with line feed (handy all along)
crlf = char(10); % [char(13),char(10)];

% % Make range vector strings (float->symbolic->fortran string)
%       Convert to symbolic format
HK_RANGE    = sym(hk_range); 
RT_RANGE    = sym(rt_range);
TREH_RANGE  = sym(treh_range);
% Generate fortran compatible strings for ranges (concatenate a bit later!)
s1 = fortran(HK_RANGE);
s2 = fortran(RT_RANGE);
s3 = fortran(TREH_RANGE);

% % Make Range vector Header Strings
% Make header declaration strings for ranges
s1_h = ['      REAL, DIMENSION(1,' num2str(length(HK_RANGE  )) ') :: HK_RANGE'];
s2_h = ['      REAL, DIMENSION(1,' num2str(length(RT_RANGE  )) ') :: RT_RANGE'];
s3_h = ['      REAL, DIMENSION(1,' num2str(length(TREH_RANGE)) ') :: TREH_RANGE'];

% % Make CEC table strings (float->symbolic->fortran string)
CEC_TABLE = sym(cec_table);  % Convert cei_table to symbolic format
s4 = fortran(CEC_TABLE);     % Export from symbolic format to Fortran code (file output)

% % Make CEC grid header declaration
s4_h = ['      REAL, DIMENSION(' num2str(length(HK_RANGE  )) , ',' ...
                                 num2str(length(RT_RANGE  )) , ',' ...
                                 num2str(length(TREH_RANGE)) ,     ...
                                 ') :: CEC_TABLE'];
                               
% % Make range size integer strings
s0   = ['      N_HK   = ' , num2str(length(HK_RANGE  )) , crlf, ...
        '      N_RT   = ' , num2str(length(RT_RANGE  )) , crlf, ...
        '      N_TREH = ' , num2str(length(TREH_RANGE))];
% % Make range size integer header declaration
s0_h = ['      INTEGER N_HK, N_RT, N_TREH'];

% Concatenate strings for all data definition!
s_arrays = [s0  , crlf, crlf ...   (size hints)
            s1  , crlf, crlf ...     (hk range)
            s2  , crlf, crlf ...     (rt range)
            s3  , crlf, crlf ...   (tppt range)
            s4  , crlf, crlf]; %     (cei grid)

% Now play the dirty trick, change array contents to replace D's with E's,
% so that constant declaration does not impose double precision!
% (this will only be ok if variable names do not contain D's!)
s_arrays(s_arrays=='D') = 'E';
        
% Concatenate strings for all declaration headers!
s_headers = [s0_h, crlf, ...       (size hints)
             s1_h, crlf, ...         (hk range)
             s2_h, crlf, ...         (rt range)
             s3_h, crlf, ...       (tppt range)
             s4_h, crlf, crlf]; %    (cei grid)


fid = fopen('4_fortan_output/CEC_DATA_precheck.f90','w');
if fid>=0
    fprintf(fid, '%s\n', s_headers);
    fprintf(fid, '%s\n', s_arrays);
    fclose(fid);
end
% % Export done!


% % Aggregate into ND-grids
cei_function = @(hk, rt, t_theta_p) interp3(hk_grid, rt_grid, treh_grid, cec_grid, hk, rt, t_theta_p, 'linear');

% % Visualize in Cuts
figure(3);
plot(hk_range , cei_function(hk_range, min(rt_range), 50), hk_range , cei_function(hk_range, mean(rt_range), 50) , hk_range , cei_function(hk_range, max(rt_range), 50))
hold on; grid on; 
xlabel('hk - Shape Factor') ; ylabel('(C_f^{rough}-C_f^{0})/C_f^{0} - Skin friction ratio (RCFR)');
title('Influence of Momentum Thickness Reynolds @ Re_h = 50')
legend(['Re_{theta} = ', num2str(min(rt_range)) ], ['Re_{theta} = ' , num2str(mean(rt_range)) ], ['Re_{theta} = ' , num2str(max(rt_range))])
% Print
set(gcf, 'PaperType', 'A5'); orient landscape
print -dpdf ./figures_CEC_export_cut_Re_h50.pdf

figure(4);
plot(hk_range , cei_function(hk_range, min(rt_range),100), hk_range , cei_function(hk_range, mean(rt_range),100) , hk_range , cei_function(hk_range, max(rt_range),100))
hold on; grid on; 
xlabel('hk - Shape Factor') ; ylabel('(C_f^{rough}-C_f^{0})/C_f^{0} - Skin friction ratio (RCFR)');
title('Influence of Momentum Thickness Reynolds @ Re_h = 100')
legend(['Re_{theta} = ', num2str(min(rt_range)) ], ['Re_{theta} = ' , num2str(mean(rt_range)) ], ['Re_{theta} = ' , num2str(max(rt_range))])
% Print
set(gcf, 'PaperType', 'A5'); orient landscape
print -dpdf ./figures_CEC_export_cut_Re_h100.pdf


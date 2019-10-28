x_over_h_list = [5 6 7 8 9 10 25 50];

% Compute integral boundary layer parameters on experimental data
% First extrapolate fields and make that default



% Allocate
rt_list         = zeros(size(x_over_h_list));
hk_list         = zeros(size(x_over_h_list));
tet_list      = zeros(size(x_over_h_list));


n_x_over_h      = 8;

target_x = x_over_h_list(n_x_over_h); 
target_z = linspace(-1,1);          % List of z/S stances
Hk_z     = zeros(size(target_z ));  % Shape factor at z/S stance
theta_z  = zeros(size(target_z ));  % Momentum thicknes at z/S stance (over delta_ref)
for n_z = 1:length(target_z)
    [~, theta_z(n_z), Hk_z(n_z)] = ED.integral_parameters(target_x,  target_z(n_z), 2.0, u_inf);
end
rt_z     = theta_z*u_inf/nu_inf;    % Reynolds theta at z/S stance

% Take reference values as averages
tet_list(n_x_over_h) = mean(theta_z);
rt_list(n_x_over_h)  = mean(rt_z   );
hk_list(n_x_over_h)  = mean(Hk_z   );

disp(['tet = ' , num2str(tet_list)])
disp(['rt  = ' , num2str(rt_list )])
disp(['hk  = ' , num2str(hk_list )])



% 
figure(1)
n_case = 1; 
subplot(221);
surf(ED.data_struct_cell{n_case}.z/ED.S, ED.data_struct_cell{n_case}.y/ED.delta_ref, ED.data_struct_cell{n_case}.rst_xx); view(2); shading flat
title(['RST_XX  -  x/h = ' , num2str(ED.d_list(n_case))]);
subplot(222);
surf(ED.data_struct_cell{n_case}.z/ED.S, ED.data_struct_cell{n_case}.y/ED.delta_ref, ED.data_struct_cell{n_case}.rst_yy); view(2); shading flat
title(['RST_YY  -  x/h = ' , num2str(ED.d_list(n_case))]);
subplot(223);
surf(ED.data_struct_cell{n_case}.z/ED.S, ED.data_struct_cell{n_case}.y/ED.delta_ref, ED.data_struct_cell{n_case}.rst_zz); view(2); shading flat
title(['RST_ZZ  -  x/h = ' , num2str(ED.d_list(n_case))]);
subplot(224);
surf(ED.data_struct_cell{n_case}.z/ED.S, ED.data_struct_cell{n_case}.y/ED.delta_ref, ED.data_struct_cell{n_case}.rst_xx + ED.data_struct_cell{n_case}.rst_yy + ED.data_struct_cell{n_case}.rst_zz); view(2); shading flat
title(['trace(RST) -  x/h = ' , num2str(ED.d_list(n_case))]);

figure(2)
n_case = 1; 
subplot(221);
surf(ED.data_struct_cell{n_case}.z/ED.S, ED.data_struct_cell{n_case}.y/ED.delta_ref, ED.data_struct_cell{n_case}.rst_xx); view(2); shading flat
title(['RST_XX  -  x/h = ' , num2str(ED.d_list(n_case))]);
subplot(222);
surf(ED.data_struct_cell{n_case}.z/ED.S, ED.data_struct_cell{n_case}.y/ED.delta_ref, ED.data_struct_cell{n_case}.rst_yy); view(2); shading flat
title(['RST_YY  -  x/h = ' , num2str(ED.d_list(n_case))]);
subplot(223);
surf(ED.data_struct_cell{n_case}.z/ED.S, ED.data_struct_cell{n_case}.y/ED.delta_ref, ED.data_struct_cell{n_case}.rst_zz); view(2); shading flat
title(['RST_ZZ  -  x/h = ' , num2str(ED.d_list(n_case))]);
subplot(224);
surf(ED.data_struct_cell{n_case}.z/ED.S, ED.data_struct_cell{n_case}.y/ED.delta_ref, ED.data_struct_cell{n_case}.rst_xx + ED.data_struct_cell{n_case}.rst_yy + ED.data_struct_cell{n_case}.rst_zz); view(2); shading flat
title(['trace(RST) -  x/h = ' , num2str(ED.d_list(n_case))]);




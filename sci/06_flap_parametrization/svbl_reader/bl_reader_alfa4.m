% BL data reader

% Set file to read
filename_root       = 'du25plr16';

% Read it
bldata = read_svbl_files(filename_root);

% % Develop x-derivative of pressure coefficient
xun_top_range   = bldata.xun_top_range;
xun_bot_range   = bldata.xun_bot_range;
cpx_top         = bldata.cpx_fun_top(5, xun_top_range);
cpx_bot         = bldata.cpx_fun_bot(5, xun_bot_range);
dcpx_dx_top     = gradient(    cpx_top, xun_top_range);
dcpx_dx_bot     = gradient(    cpx_bot, xun_bot_range);

sqrt(trapz(xun_top_range, dcpx_dx_top.^2)) + sqrt(trapz(xun_bot_range, dcpx_dx_bot.^2))

figure(1); plot(xun_top_range, cpx_top    ); grid on; hold on;
figure(1); plot(xun_bot_range, cpx_bot    ); grid on; hold on;
figure(2); plot(xun_top_range, dcpx_dx_top); grid on;  hold on;
figure(2); plot(xun_bot_range,-dcpx_dx_bot); grid on;  hold on;


xun_top_range   = bldata.xun_top_range;
xun_bot_range   = bldata.xun_bot_range;
cpx_top_rfun    = @(aoa)      bldata.cpx_fun_top(aoa, xun_top_range);
cpx_bot_rfun    = @(aoa)      bldata.cpx_fun_bot(aoa, xun_bot_range);
dcpx_dx_top_rfun= @(aoa)      gradient(cpx_top_rfun(aoa), xun_top_range);
dcpx_dx_bot_rfun= @(aoa)      gradient(cpx_bot_rfun(aoa), xun_bot_range);
dcpx_dx_top_fun = @(aoa, xun) interp1(xun_top_range, dcpx_dx_top_rfun(aoa), xun);
dcpx_dx_bot_fun = @(aoa, xun) interp1(xun_bot_range, dcpx_dx_bot_rfun(aoa), xun);

RMS_dcpx_dx_fun = @(aoa)      sqrt(trapz(xun_top_range, dcpx_dx_top_rfun(aoa).^2)) + sqrt(trapz(xun_bot_range, dcpx_dx_bot_rfun(aoa).^2));

% % Some Plotting
figure(11)
plot(bldata.aoa_range, bldata.cpx_fun_top(bldata.aoa_range, 0.05)); grid on; hold on;
plot(bldata.aoa_range, bldata.cpx_fun_top(bldata.aoa_range, 0.20)); grid on;
plot(bldata.aoa_range, bldata.cpx_fun_top(bldata.aoa_range, 0.40)); grid on;
plot(bldata.aoa_range, bldata.cpx_fun_top(bldata.aoa_range, 0.60)); grid on;
plot(bldata.aoa_range, bldata.cpx_fun_top(bldata.aoa_range, 0.80)); grid on;
plot(bldata.aoa_range, bldata.cpx_fun_top(bldata.aoa_range, 0.95)); grid on;
xlabel('\alpha - Angle of attack'); ylabel('C_p - Pressure Coefficient');
title('Pressure Coefficient vs Angle of Attack at Chordwise Stances (Top Side)');
legend(['x/c = ' num2str(0.05)] , ...
       ['x/c = ' num2str(0.20)] , ...
       ['x/c = ' num2str(0.40)] , ...
       ['x/c = ' num2str(0.60)] , ...
       ['x/c = ' num2str(0.80)] , ...
       ['x/c = ' num2str(0.95)]);
   
figure(12)
plot(bldata.aoa_range, bldata.h_fun_top(bldata.aoa_range, 0.05)); grid on; hold on;
plot(bldata.aoa_range, bldata.h_fun_top(bldata.aoa_range, 0.20)); grid on;
plot(bldata.aoa_range, bldata.h_fun_top(bldata.aoa_range, 0.40)); grid on;
plot(bldata.aoa_range, bldata.h_fun_top(bldata.aoa_range, 0.60)); grid on;
plot(bldata.aoa_range, bldata.h_fun_top(bldata.aoa_range, 0.80)); grid on;
plot(bldata.aoa_range, bldata.h_fun_top(bldata.aoa_range, 0.95)); grid on;
xlabel('\alpha - Angle of attack'); ylabel('H - Shape Factor');
title('Shape Factor vs Angle of Attack at Chordwise Stances');
legend(['x/c = ' num2str(0.05)] , ...
       ['x/c = ' num2str(0.20)] , ...
       ['x/c = ' num2str(0.40)] , ...
       ['x/c = ' num2str(0.60)] , ...
       ['x/c = ' num2str(0.80)] , ...
       ['x/c = ' num2str(0.95)]);

   
figure(13)
plot(bldata.aoa_range, bldata.cfx_fun_top(bldata.aoa_range, 0.05)); grid on; hold on;
plot(bldata.aoa_range, bldata.cfx_fun_top(bldata.aoa_range, 0.20)); grid on;
plot(bldata.aoa_range, bldata.cfx_fun_top(bldata.aoa_range, 0.40)); grid on;
plot(bldata.aoa_range, bldata.cfx_fun_top(bldata.aoa_range, 0.60)); grid on;
plot(bldata.aoa_range, bldata.cfx_fun_top(bldata.aoa_range, 0.80)); grid on;
plot(bldata.aoa_range, bldata.cfx_fun_top(bldata.aoa_range, 0.95)); grid on;
xlabel('\alpha - Angle of attack'); ylabel('C_f - Skin Friction Coefficient');
title('Skin Friction Coefficient vs Angle of Attack at Chordwise Stances');
legend(['x/c = ' num2str(0.05)] , ...
       ['x/c = ' num2str(0.20)] , ...
       ['x/c = ' num2str(0.40)] , ...
       ['x/c = ' num2str(0.60)] , ...
       ['x/c = ' num2str(0.80)] , ...
       ['x/c = ' num2str(0.95)]);

   

% % Some plotting (Cp)
% aoa = 8;
% xun_plot = linspace(1, 0);
% plot(xun_plot, -base_cpx_fun(aoa, 1-xun_plot));
% hold on; grid on;
% plot(  xun_plot, -base_cpx_fun(aoa, xun_plot+1));
% legend('Upper Side', 'Lower Side');
% % Some plotting (Cf)
% aoa = 8;
% plot(xun_plot, base_cfx_fun(aoa, 1-xun_plot));
% hold on; grid on;
% plot(  xun_plot, base_cfx_fun(aoa, xun_plot+1));
% legend('Upper Side', 'Lower Side');
% 
% 
% % New Plotting!
% aoa = 8;
% xun_plot = linspace(1, 0, 1000);
% plot(xun_plot, cpx_fun_top(aoa, xun_plot));
% hold on; grid on;
% plot(  xun_plot, cpx_fun_bot(aoa, xun_plot));
% legend('Upper Side', 'Lower Side');
% % Validation
% plot(base_bldata_xun, -unique_sorted_bldata_cpx( :, 61), 'x')



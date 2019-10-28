        % % Plot Velocity Profiles
        % Centerline
        set(0, 'CurrentFigure', 14)
        subplot(231)
        j_center = CM.j_center;
        h4_231   = plot( new_pod_state.u_tilde(:, j_center )/ new_pod_state.ue, CM.y_range/delta0, ...
                                   u_bar(:, j_center )/ new_pod_state.ue, CM.y_range/delta0, ...
                                   u_bar(:, j_center )/ new_pod_state.ue + new_pod_state.u_tilde(:, j_center )/ new_pod_state.ue, CM.y_range/delta0, ...
                                   ED.interp_velocities(state.x / h_VG, CM.y_range/ delta0, 0) / u_inf, CM.y_range/delta0, 'k');
        grid on; axis([-0.3 1.1 0 1.5]); %grid on; axis([-0.5 1.5 0 2]);
        % title('Central Symmetry Line');
        title(['X/h = ' num2str(round(state.x / h_VG*100)/100), '     Z/S = ', num2str(0)])
        xlabel('U/U_{e}'); ylabel('Y/\delta_{ref}');
        % Mid Line
        subplot(233)
        j_quarter = round(j_center/2);
        h4_233    = plot(new_pod_state.u_tilde(:, j_quarter)/ new_pod_state.ue, CM.y_range/delta0, ...
                                   u_bar(:, j_quarter)/ new_pod_state.ue, CM.y_range/delta0, ...
                                   u_bar(:, j_quarter)/ new_pod_state.ue + new_pod_state.u_tilde(:, j_quarter)/ new_pod_state.ue, CM.y_range/delta0, ...
                                   ED.interp_velocities(state.x / h_VG, CM.y_range/ delta0, -1) / u_inf, CM.y_range/delta0, 'k');
        grid on; axis([-0.3 1.1 0 1.5]); %grid on; axis([-0.5 1.5 0 2]);
        % title('Side Symmetry Line');
        title(['X/h = ' num2str(round(state.x / h_VG*100)/100), '     Z/S = ', num2str(-1)])
        xlabel('U/U_{e}'); ylabel('Y/\delta_{ref}');
        % legend('Mixed Flow' , 'Pure Shear Flow', 'Total Predicted Flow', 'Experimental Result')
        legend(h4_233([4 3 2 1]), 'Experiment' , 'Prediction', 'Shear Component', 'Mixing Component')
        % Intermediate Line
        subplot(232)
        j_eighth = round((j_center+j_quarter)/2);
        h4_232   = plot(new_pod_state.u_tilde(:, j_eighth)/ new_pod_state.ue, CM.y_range/delta0, ...
                                  u_bar(:, j_eighth)/ new_pod_state.ue, CM.y_range/delta0, ...
                                  u_bar(:, j_eighth)/ new_pod_state.ue + new_pod_state.u_tilde(:, j_eighth)/ new_pod_state.ue, CM.y_range/delta0, ...
                                  ED.interp_velocities(state.x / h_VG, CM.y_range/ delta0, - CM.z_range(j_eighth) ./ CM.S) / u_inf, CM.y_range/delta0, 'k');
        grid on; axis([-0.3 1.1 0 1.5]); % grid on; axis([-0.5 1.5 0 2]);
        % title('Intermediate Line');
        title(['X/h = ' num2str(round(state.x / h_VG*100)/100), '     Z/S = ', num2str(0.5)])
        xlabel('U/U_{e}'); ylabel('Y/\delta_{ref}');
        % Update colors
        blue   = [0         0.4470    0.7410];
        red    = [0.8500    0.3250    0.0980];
        yellow = [0.9290    0.6940    0.1250];
        black  = [0         0         0     ];
        h4_231(1).Color = black ; h4_232(1).Color = h4_231(1).Color; h4_233(1).Color = h4_231(1).Color;
        h4_231(2).Color = yellow; h4_232(2).Color = h4_231(2).Color; h4_233(2).Color = h4_231(2).Color;
        h4_231(3).Color = red   ; h4_232(3).Color = h4_231(3).Color; h4_233(3).Color = h4_231(3).Color;
        h4_231(4).Color = blue  ; h4_232(4).Color = h4_231(4).Color; h4_233(4).Color = h4_231(4).Color;
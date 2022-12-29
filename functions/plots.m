
%% De-tumbling
if alg_index == 1 || alg_index == 3

    % De-tumbling angular velocities
      phase1.tout = phase1.tout/3600; % convert time from seconds to hours
    abs_w = zeros(length(phase1.tout),1);
    for i = 1:length(phase1.tout)
        abs_w(i) = norm(phase1.w_BN(i,:));
    end

    fig_phase1_w = figure('Units','normalized','OuterPosition',[0 0 1 1]);
    hold on;
    grid on;
    plot(phase1.tout,phase1.w_BN(:,1))
    plot(phase1.tout,phase1.w_BN(:,2))
    plot(phase1.tout,phase1.w_BN(:,3))
    plot(phase1.tout,abs_w)

    xline(sensors.star.startup/3600,'--','LineWidth',1.5);
    text(sensors.star.startup/3600+0.1,2.4,'Star sensors minimum start up time','FontSize',12)

    xline(phase1.tout(end),'--','LineWidth',1.5);
    text(phase1.tout(end)-0.75,2.4,'Angular velocity requirements satisfaction time','FontSize',12)

    legend('$\omega_x$','$\omega_y$','$\omega_z$','$||\omega||$','','','Location','northeast');

    xlabel('Elapsed time [h]')
    ylabel('Angular velocity [deg/s]')

    % De-tumbling Control Torques
    fig_phase1_T = figure('Units','normalized','OuterPosition',[0 0 1 1]);
    hold on;
    grid on;
    plot(phase1.tout,phase1.applied_torque(:,1))
    plot(phase1.tout,phase1.applied_torque(:,2))
    plot(phase1.tout,phase1.applied_torque(:,3))

    legend('$T_x$','$T_y$','$T_z$','Location','northeast');

    xlabel('Elapsed time [h]')
    ylabel('Applied control torque [Nm]')

    % save plots
    if strcmp(save_plots,'yes')
        %de-tumbling
        exportgraphics(fig_phase1_w,'.\report_images\w_detumbling.eps','ContentType','vector')
        exportgraphics(fig_phase1_T,'.\report_images\T_detumbling.eps','ContentType','vector')
    end

end

%% Pointing
if alg_index == 2 || alg_index == 3

    %Pointing angular velocities
    phase2.tout = phase2.tout/3600; % convert time from seconds to hours

    fig_phase2_w = figure('Units','normalized','OuterPosition',[0 0 1 1]);
    hold on;
    grid on;
    plot(phase2.tout,phase2.w_BN(:,1))
    plot(phase2.tout,phase2.w_BN(:,2))
    plot(phase2.tout,phase2.w_BN(:,3))

    legend('$\omega_x$','$\omega_y$','$\omega_z$','Location','northeast');

    xlabel('Elapsed time [h]')
    ylabel('Angular velocity [deg/s]')

    %Pointing Control Torques
    fig_phase2_T = figure('Units','normalized','OuterPosition',[0 0 1 1]);
    hold on;
    grid on;
    plot(phase2.tout,phase2.applied_torque(:,1))
    plot(phase2.tout,phase2.applied_torque(:,2))
    plot(phase2.tout,phase2.applied_torque(:,3))

    legend('$T_x$','$T_y$','$T_z$','Location','northeast');

    xlabel('Elapsed time [h]')
    ylabel('Applied control torque [Nm]')

    % Pointing Euler angles
    fig_phase2_E = figure('Units','normalized','OuterPosition',[0 0 1 1]);
    hold on;
    grid on;
    plot(phase2.tout,phase2.E_312(:,1))
    plot(phase2.tout,phase2.E_312(:,2))
    plot(phase2.tout,phase2.E_312(:,3))

    legend('$\phi$','$\vartheta$','$\psi$','Location','northeast');

    xlabel('Elapsed time [h]')
    ylabel('Euler angles [deg]')

    % Pointing star sensors
    phase2.sensors_status.Time = phase2.sensors_status.Time/3600;
    phase2.sensors_status.Data = squeeze(phase2.sensors_status.Data);
    fig_phase2_sens = figure('Units','normalized','OuterPosition',[0 0 1 1]);
    hold on;
    grid on;
    plot(phase2.sensors_status.Time,phase2.sensors_status.Data)

    xlabel('Elapsed time [h]')
    ylabel('status [1=on 0=off]')

    % Pointing error
    phase2.Pointing_error.Time = phase2.Pointing_error.Time/3600;
    fig_phase2_err = figure('Units','normalized','OuterPosition',[0 0 1 1]);
    hold on;
    grid on;
    plot(phase2.Pointing_error.Time,phase2.Pointing_error.Data(:,1))
    plot(phase2.Pointing_error.Time,phase2.Pointing_error.Data(:,2))
    plot(phase2.Pointing_error.Time,phase2.Pointing_error.Data(:,3))

    if phase2.tout(end)<Simtime
        xline(phase2.Pointing_error.Time(end)-settings.Time/3600,'--','LineWidth',1.5);
        text(phase2.Pointing_error.Time(end)-settings.Time/3600+0.2,0.55,'End of slew maneouvre','FontSize',12)
    end
    legend('$A_e^{31}$','$A_e^{13}$','$A_e^{21}$','','Location','northeast')
    xlabel('Elapsed time [h]')
    ylabel('$A_e$')

    % save plots
    if strcmp(save_plots,'yes')

        exportgraphics(fig_phase2_w,'.\report_images\w_pointing.eps','ContentType','vector')
        exportgraphics(fig_phase2_T,'.\report_images\T_pointing.eps','ContentType','vector')
        exportgraphics(fig_phase2_E,'.\report_images\E_pointing.eps','ContentType','vector')
        exportgraphics(fig_phase2_sens,'.\report_images\sensors_pointing.eps','ContentType','vector')
        exportgraphics(fig_phase2_err,'.\report_images\pointing.eps','ContentType','vector')
    end

end

%% No control 
if alg_index == 4
    
 % No control Euler angles
 phase0.tout = phase0.tout/3600;
    fig_phase0_E = figure('Units','normalized','OuterPosition',[0 0 1 1]);
    hold on;
    grid on;
    plot(phase0.tout,phase0.E_312(:,1))
    plot(phase0.tout,phase0.E_312(:,2))
    plot(phase0.tout,phase0.E_312(:,3))

    legend('$\phi$','$\vartheta$','$\psi$','Location','northeast');

    xlabel('Elapsed time [h]')
    ylabel('Euler angles [deg]')

    %angular velocities
    fig_phase0_w = figure('Units','normalized','OuterPosition',[0 0 1 1]);
    hold on;
    grid on;
    plot(phase0.tout,phase0.w_BN(:,1))
    plot(phase0.tout,phase0.w_BN(:,2))
    plot(phase0.tout,phase0.w_BN(:,3))

    legend('$\omega_x$','$\omega_y$','$\omega_z$','Location','northeast');

    xlabel('Elapsed time [h]')
    ylabel('Angular velocity [deg/s]')

      % save plots
    if strcmp(save_plots,'yes')
        exportgraphics(fig_phase0_w,'.\report_images\w_nocontrol.eps','ContentType','vector')
        exportgraphics(fig_phase0_E,'.\report_images\E_nocontrol.eps','ContentType','vector')
    end


end



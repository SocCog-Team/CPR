addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/PSY_Setup/random_walk/
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/

close all

% Parameters
N = 7200;               % Number of time steps (7200 frames, 1min)
dt = 1/120;             % Time step (120Hz)
tau_sec = .5;          % Time constant for the angular velocity decay
sigma_rpm = 2;        % Standard deviation for the random noise in omega
omega = zeros(1, N);    % Array to store angular velocities
angle = zeros(1, N);    % Array to store angles

% tau: A larger value keeps the angular velocity persistent longer, while a smaller value makes it decay faster.
% sigma: A larger value increases the randomness in the angular velocity.
% dt: Controls the resolution of the time steps.

% Initial conditions
omega(1) = 2 * pi * (rand - 0.5);  % Initial random angular velocity
angle(1) = 2 * pi * rand;          % Initial random angle between 0 and 2*pi

% Generate the random walk
for i = 2:N
    % Gaussian random noise with standard deviation scaled by sqrt(dt)
    GRN = randn;
    
    % Update omega [Revolutions per minute]
    omega(i) = omega(i-1) - (omega(i-1) / tau_sec) * dt + sigma_rpm * GRN * sqrt(2*dt/tau_sec);
    
    % Update angle(t+dt)
    angle(i) = angle(i-1) + omega(i) * dt;
    
    % Keep the angle within [0, 2*pi]
    angle(i) = mod(angle(i), 2 * pi);
end

figure;
polarplot(angle, 1:N);
title('Random Walk with Inertia in Polar Coordinates');

% Write 2 File 4 MWorks
out                         = CPR_create_random_walk_v2();
out.RDP_direction           = rad2deg(angle);
[~]                     	= CPR_write_txt(out,'/Users/fschneider/Desktop/');  

%% PLOT
%%% Direction and events
f = figure('units','pixels','position',[0 300 560 420]); hold on
plot(out.RDP_direction, 'k','LineWidth', 1.5) % Plot direction
scatter(out.feedback_ts, out.RDP_direction(out.feedback_ts), 'ro', 'filled') % Plot reward
% scatter(out.jump_ts, rad2deg(out.dir_seed), 'b^', 'filled') % Plot direction of random seeds

xlim([0 7200])
ylim([0 360])
xlabel('Frames')
ylabel('Direction [deg]')
set(gca, 'fontsize', 14)
% legend('Dir','Rew','Seed')
legend('Dir','Rew')
  
%%% Polar histogram
f = figure('units','pixels','position',[562 300 560 420]); hold on; clf
ax                          = polaraxes;
pp                          = polarhistogram(deg2rad(out.RDP_direction), 180);
pp.FaceColor                = [.2 .2 .2];
pp.EdgeColor                = 'none';
ax.ThetaZeroLocation        = 'top';
ax.ThetaDir                 = 'clockwise';
ax.FontSize                 = 14;

%%% Circular difference
f = figure('units','pixels','position',[953 300 560 420]); hold on
pp                          = plot(rad2deg(circ_dist(deg2rad(out.RDP_direction(1:end-1)),deg2rad(out.RDP_direction(2:end)))));
ax                          = gca;
ax.YLabel.String            = 'Difference (deg)';
ax.XLabel.String            = 'Samples [#]';
ax.FontSize                 = 14;
ax.XLim                     = [0 7200];

%%
f = figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
for iSim = 1:16
    % Create new stimulus
    out = CPR_create_random_walk_v2();
    dirs{iSim} = out.RDP_direction;
    
    % Visualise
    subplot(4,4,iSim); hold on
    plot(out.RDP_direction, 'k','LineWidth', 1.5) % Plot direction
    scatter(out.feedback_ts, out.RDP_direction(out.feedback_ts), 'ro', 'filled') % Plot reward
%     scatter(out.jump_ts, rad2deg(out.dir_seed), 'b^', 'filled') % Plot direction of random seeds
    
    xlim([0 7200])
    ylim([0 360])
    xlabel('Frames')
    ylabel('Direction [deg]')
    set(gca, 'fontsize', 14)
    set(gca, 'xtick', [0 3600 7200])
    set(gca, 'ytick', [0 180 360])
end
% legend('Dir','Rew','Seed')
legend('Dir','Rew')
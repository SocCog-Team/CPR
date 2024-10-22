addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/PSY_Setup/random_walk/
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/

close all

%%% Generate new stimulus
out = CPR_create_random_walk_v2();

%%% Direction and events
f = figure('units','pixels','position',[0 300 560 420]); hold on
plot(out.RDP_direction, 'k','LineWidth', 1.5) % Plot direction
scatter(out.feedback_ts, out.RDP_direction(out.feedback_ts), 'ro', 'filled') % Plot reward
scatter(out.jump_ts, rad2deg(out.dir_seed), 'b^', 'filled') % Plot direction of random seeds

xlim([0 7200])
ylim([0 360])
xlabel('Frames')
ylabel('Direction [deg]')
set(gca, 'fontsize', 14)
legend('Dir','Rew','Seed')
  
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
% f = figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]);
% for iSim = 1:16
%     % Create new stimulus
%     out = CPR_create_random_walk_v2();
%     dirs{iSim} = out.RDP_direction;
%     
%     % Visualise
%     subplot(4,4,iSim); hold on
%     plot(out.RDP_direction, 'k','LineWidth', 1.5) % Plot direction
%     scatter(out.feedback_ts, out.RDP_direction(out.feedback_ts), 'ro', 'filled') % Plot reward
%     scatter(out.jump_ts, rad2deg(out.dir_seed), 'b^', 'filled') % Plot direction of random seeds
%     
%     xlim([0 7200])
%     ylim([0 360])
%     xlabel('Frames')
%     ylabel('Direction [deg]')
%     set(gca, 'fontsize', 14)
%     set(gca, 'xtick', [0 3600 7200])
%     set(gca, 'ytick', [0 180 360])
% end
% legend('Dir','Rew','Seed')
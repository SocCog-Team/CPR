addpath /Users/fschneider/Documents/MATLAB/CircStat2012a
clear all
close all

a = [350 0 90 180 270 360 90 180 90 0 270 180]; % deg
b = deg2rad(a); % rad

for i = 1:length(a)-1
    dff(i) = rad2deg(circ_dist(b(i),b(i+1))); % Difference [deg]
end

figure
ax1 = subplot(2,1,1);
stairs(a)
ax1.XLim = [1 length(a)];
ax1.YLabel.String = 'Direction [deg]';
ax1.FontSize = 16;

ax2 = subplot(2,1,2);
stairs(cumsum(-dff))
ax2.XLim = [1 length(a)];
ax2.YLabel.String = 'Difference [deg]';
ax2.FontSize = 16;


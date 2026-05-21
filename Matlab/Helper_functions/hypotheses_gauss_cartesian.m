x = linspace(-5, 5, 1000);

% Base gaussian
sigma = 1; mu = 0;
g = @(x, mu, s, a) a * exp(-((x-mu).^2)/(2*s^2));

figure;

% --- Left panel: Multiplicative scaling ---
ax1 = subplot(1,3,1);
plot(x, g(x,0,1,0.4), 'g', 'LineWidth', 2); hold on;
plot(x, g(x,0,1,1.0), 'k', 'LineWidth', 2);
plot(x, g(x,0,1,1.6), 'r', 'LineWidth', 2);
title('Multiplicative Scaling');
xlabel('Firing rate'); ylabel('Stimulius direction');
xlim([-5 5]); ylim([0 2]);

% --- Center panel: Additive changes ---
ax2 = subplot(1,3,2);
base = g(x,0,1,1.0);
plot(x, base, 'g', 'LineWidth', 2); hold on;
plot(x, base + 0.3,        'k', 'LineWidth', 2);
plot(x, base + 0.6, 'r', 'LineWidth', 2);
title('Additive Offset');
xlabel('x');
xlim([-5 5]); ylim([0 2]);

% --- Right panel: Width changes ---
ax3 = subplot(1,3,3);
plot(x, g(x,0,0.5,1), 'g', 'LineWidth', 2); hold on;
plot(x, g(x,0,1.0,1), 'k', 'LineWidth', 2);
plot(x, g(x,0,1.8,1), 'r', 'LineWidth', 2);
title('Width Changes');
xlabel('x');
xlim([-5 5]); ylim([0 2]);

% Make all subpanels square
for ax = [ax1 ax2 ax3]
    axis(ax, 'square');
end

box off

print(gcf,'/Users/fschneider/Desktop/gauss', '-r300', '-dsvg');shg
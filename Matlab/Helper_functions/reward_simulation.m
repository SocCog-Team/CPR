addpath /Users/fschneider/Documents/MATLAB/cbrewer/

clear all
close all

rew_power       = [.2:.2:2];
js_tlt          = [0:.01:1];
rew_max_ml      = .5;
rew_min_ml      = .05;
[col]           = cbrewer('div', 'RdYlBu', length(rew_power));

mhir            = [ 0.2153    0.5317    0.6465    0.6939]; % based on data
coh             = [0 .39 .56 .98]; % based on paradigm
p_rew           = polyfit(coh, mhir, 2);
y_fit           = polyval(p_rew, linspace(0, 1, 101));


f1 = figure;hold on
f2 = figure;hold on
for i = 1:length(rew_power)
    rew_target_ml   = rew_min_ml + ((js_tlt.^rew_power(i)) * (rew_max_ml - rew_min_ml));
    
    figure(f1)
    p_rew(i)        = plot(rew_target_ml,'LineWidth', 2, 'Color',col(i,:));
    str{i}          = num2str(rew_power(i));

    figure(f2)
    p_hit(i)        = plot(rew_target_ml .* y_fit,'LineWidth', 2, 'Color',col(i,:));
end

figure(f1)
xlim([0 100])
xlabel('joystick tilt')
ylabel('reward volume [ml] | hit')
title({'Payoff functions based on exponents';'Rew [ml] = minRew + ((jsTILT ^^ rewPOWER) * (maxREW-minRew))'})
legend(p_rew,str,'Location','southeast')

figure(f2)
xlim([0 100])
xlabel('coherence == tilt')
ylabel('avg. reward volume [ml] / hit')
title({'Payoff scaled by hit rate (assumes tilt == coh)'; 'reward function .* hit rate'})
legend(p_hit,str,'Location','southeast')


% Within coherence
% modulation of tilt affect hit rate
% What about the overall level of reward
% For each reward function

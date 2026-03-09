close all
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab

%%
c = 0;
for iSubj = 1:length(solo_perf)

    if isempty(dyad_perf{iSubj})
        continue
    end
    c = c+1;
    acc_dff(c) = mean(dyad_perf{iSubj}.macc_state - solo_perf{iSubj}.macc_state);
    tlt_dff(c) = mean(dyad_perf{iSubj}.mecc_state - solo_perf{iSubj}.mecc_state);
end

f               = figure; hold on
vl              = violinplot([tlt_dff',acc_dff']);
ln              = line([0 3], [0 0],'LineStyle',':','LineWidth',2, 'Color','k');

ax              = gca;
ax.FontSize     = 20;
ax.XTickLabel   = {'Confidence', 'Accuracy'};
ax.YLabel.String= 'Difference (Dyad - Solo)';
[p_tlt,h]       = signrank(tlt_dff);
[p_acc,h]       = signrank(acc_dff);

text(.7, -.05 , ['p = ' num2str(p_tlt)], 'FontSize', 14, 'Color', 'k')
text(1.7, -.05 , ['p = ' num2str(p_acc)], 'FontSize', 14, 'Color', 'k')
text(.7, -.06 , ['median = ' num2str(median(tlt_dff))], 'FontSize', 14, 'Color', 'k')
text(1.7, -.06 , ['median = ' num2str(median(acc_dff))], 'FontSize', 14, 'Color', 'k')
print(gcf, ['/Users/fschneider/Desktop/raw_dff_avg'], '-r500', '-dsvg', '-painters');

%% Raw difference: dyad minus solo
tlt            = [raw.decc1 - raw.ecc1 raw.decc2 - raw.ecc2]; % joystick tilt
acc            = [raw.dacc1 - raw.acc1 raw.dacc2 - raw.acc2]; % joystick accuracy

f               = figure; hold on
vl              = violinplot([tlt',acc']);
ln              = line([0 3], [0 0],'LineStyle',':','LineWidth',2, 'Color','k');

ax              = gca;
ax.FontSize     = 20;
ax.XTickLabel   = {'Confidence', 'Accuracy'};
ax.YLabel.String= 'Difference (Dyad - Solo)';

[p_tlt,h]       = signrank(tlt);
[p_acc,h]       = signrank(acc);

text(.7, -.18 , ['p = ' num2str(p_tlt)], 'FontSize', 14, 'Color', 'k')
text(1.7, -.18 , ['p = ' num2str(p_acc)], 'FontSize', 14, 'Color', 'k')
text(.7, -.16 , ['median = ' num2str(median(tlt))], 'FontSize', 14, 'Color', 'k')
text(1.7, -.16 , ['median = ' num2str(median(acc))], 'FontSize', 14, 'Color', 'k')
print(gcf, ['/Users/fschneider/Desktop/raw'], '-r500', '-dsvg', '-painters');

%% Convergence

%%% Tilt
p1_better               = raw.ecc1 > raw.ecc2;
ecc_solo_dff            = [raw.ecc1(p1_better) - raw.ecc2(p1_better) raw.ecc2(~p1_better) - raw.ecc1(~p1_better)];
ecc_dyad_dff            = [raw.decc1(p1_better) - raw.decc2(p1_better) raw.decc2(~p1_better) - raw.decc1(~p1_better)];
plot_this_data(ecc_solo_dff, ecc_dyad_dff,'joystick tilt',-.05)

%%% Accuracy
p1_better               = raw.acc1 > raw.acc2;
acc_solo_dff            = [raw.acc1(p1_better) - raw.acc2(p1_better) raw.acc2(~p1_better) - raw.acc1(~p1_better)];
acc_dyad_dff            = [raw.dacc1(p1_better) - raw.dacc2(p1_better) raw.dacc2(~p1_better) - raw.dacc1(~p1_better)];
plot_this_data(acc_solo_dff, acc_dyad_dff,'joystick accuracy',-.02)


function plot_this_data(in_solo, in_dyad,str,yofs)
col = [.6 .1 .1; .6 .9 .9];
nbin = 10;

figure
vl                  = violinplot([in_solo',in_dyad']);
vl                  = improveViolin(vl,col);
set(gca,'FontSize',20)
set(gca,'XTickLabel',{'solo', 'dyad'})
ylabel('Inter-player difference (better - worse)')
title(str)
ln = line([0 2], [0 0],'LineStyle',':')
print(gcf, ['/Users/fschneider/Desktop/js_dff_' str], '-r500', '-dsvg', '-painters');

[p,h] = signrank(in_solo,in_dyad)
text(1.25, yofs , ['p = ' num2str(p)], 'FontSize', 14, 'Color', 'k')

figure
histogram(in_solo, nbin, 'DisplayStyle', 'stairs','EdgeColor', col(1,:), 'LineWidth', 5);
hold on
histogram(in_dyad, nbin, 'DisplayStyle', 'stairs','EdgeColor', col(2,:), 'LineWidth', 5);
set(gca,'FontSize',20)
ylabel('#')
xlabel('Inter-player difference (better - worse)')
title(str)
end


function vl = improveViolin(vl,col_map)
for iV=1:length(vl)
    vl(iV).BoxWidth                     = .025;
    vl(iV).ViolinColor{1}               = col_map(iV,:);
    vl(iV).ViolinAlpha{1}               = .8;
    vl(iV).ScatterPlot.MarkerFaceColor  = [.25 .25 .25];
    vl(iV).ScatterPlot.MarkerEdgeColor  = 'none';
    vl(iV).ScatterPlot.SizeData         = 60;
end
end
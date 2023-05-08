f                   = figure; hold on
addpath /Users/fschneider/Documents/MATLAB/cbrewer
pth                 = '/Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_psychophysics/';
fname               = 'Subjects_summary.xlsx';
x                   = readtable([pth fname]);

dte                 = [x.DateTraining x.DateTraining2 x.DateSolo1 x.DateSolo2 x.DateAgent x.DateDyadic1 x.DateDyadic2 x.DateDyadic3  x.DateDyadic4];
dte(sum(isnat(dte),2) == size(dte,2),:) = [];

[dte_s, idx]        = sort(dte,2);
idx(isnat(dte_s))   = 0;
idx_tmp             = [idx zeros(size(idx,1),1)];
idx_pad             = [idx_tmp; zeros(1, size(idx_tmp,2))];
pc                  = pcolor(idx_pad);

cl_nat              = [0 0 0];
cl_train            = [.5 .5 .5];
cl_solo             = [1 0 0];
cl_agent            = [0 1 0];
cl_dyadic           = [0 0 1];

ax                  = gca;
ax.YTick            = [1:size(dte,1)] + .5;
ax.YTickLabel       = x.Abbreviation;
ax.XLabel.String    = '# Session';
ax.XTick            = [];
ax.FontSize         = 16;

col                 = cbrewer('qual', 'Dark2', 4, 'PCHIP');
cmap                = [0 0 0; col(1,:); col(1,:); col(2,:); col(2,:); col(3,:); col(4,:); col(4,:); col(4,:); col(4,:)];

% sbj_exit = [2 13];
% for iSbj = 1:length(sbj_exit)
%     id = sbj_exit(iSbj);
% plot([1 2 nan 1 2], [id id+1 nan id+1 id],'Color',[1 1 1],'LineWidth',2);
% end

h                   = zeros(5, 1);
h(1)                = plot(NaN,NaN,'.','Color', cmap(1,:));
h(2)                = plot(NaN,NaN,'.','Color', cmap(2,:));
h(3)                = plot(NaN,NaN,'.','Color', cmap(4,:));
h(4)                = plot(NaN,NaN,'.','Color', cmap(6,:));
h(5)                = plot(NaN,NaN,'.','Color', cmap(7,:));
[lg,icons]          = legend(h, 'TBD','Training/Pretest','Solo', 'Dyadic_Computer', 'Dyadic_Human', 'Interpreter', 'none');
lg.Location         = 'eastoutside';
lg.Box              = 'off';
icons               = findobj(icons,'Type','line');
icons               = findobj(icons,'Marker','none','-xor');

colormap(cmap)
set(icons,'MarkerSize',35);

axis image
axis ij

print(f, '/Users/fschneider/Desktop/exp_order', '-r300', '-dpng');

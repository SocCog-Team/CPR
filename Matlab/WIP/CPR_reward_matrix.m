close all

% Reward parameter
Ra                  = 0:.001:1;
Rb                  = 0:.001:1;
R                   = 1;

%% NEUTRAL
for iA = 1:length(Ra)
    for iB = 1:length(Rb)
        RA(iA,iB) 	= (Ra(iA) * Rb(iB)) * R;
        RB(iA,iB) 	= (Ra(iA) * Rb(iB)) * R;
    end
end

f                   = figure;
im                  = surf(Ra,Rb,RA);
im.LineStyle        = 'none';
ax                  = gca;
ax.Title.String     = 'NEUTRAL: Reward Player A';
ax.XTick            = [0 .5 1];
ax.YTick            = [0 .5 1];
ax.ZTick            = [0 .5 1];
ax.ZTickLabel       = {'0','50','100'};
ax.FontSize         = 16;
ax.XLabel.String    = 'Accuracy [norm]';
ax.YLabel.String    = 'Confidence [norm]';
ax.ZLabel.String    = 'Reward [%]';
colormap turbo
print(f, '/Users/fschneider/Desktop/neutralA', '-r300', '-dpng');

f                   = figure;
im                  = surf(Ra,Rb,RB);
im.LineStyle        = 'none';
ax                  = gca;
ax.Title.String     = 'NEUTRAL: Reward Player B';
ax.XTick            = [0 .5 1];
ax.YTick            = [0 .5 1];
ax.ZTick            = [0 .5 1];
ax.ZTickLabel       = {'0','50','100'};
ax.FontSize         = 16;
ax.XLabel.String    = 'Accuracy [norm]';
ax.YLabel.String    = 'Confidence [norm]';
ax.ZLabel.String    = 'Reward [%]';
colormap turbo
print(f, '/Users/fschneider/Desktop/neutralB', '-r300', '-dpng');

%% COMPETITION

%%% FINITE AMOUNT [Igors formula] %%%
for iA = 1:length(Ra)
    for iB = 1:length(Rb)
        RA(iA,iB) 	= (max([Ra(iA) Rb(iB)]) * R) * (Ra(iA) / (Ra(iA) + Rb(iB)));
        RB(iA,iB) 	= (max([Ra(iA) Rb(iB)]) * R) * (Rb(iB) / (Ra(iA) + Rb(iB)));
    end
end

f                   = figure;
im                  = surf(Ra,Rb,RA);
im.LineStyle        = 'none';
ax                  = gca;
ax.Title.String     = 'FINITE: Reward Player A';
ax.XTick            = [0 .5 1];
ax.YTick            = [0 .5 1];
ax.ZTick            = [0 .5 1];
ax.FontSize         = 16;
ax.XLabel.String    = {'Score B'; '[Acc * Conf]'};
ax.YLabel.String    = {'Score A'; '[Acc * Conf]'};
ax.ZLabel.String    = 'Reward [a.u.]';
colormap turbo
print(f, '/Users/fschneider/Desktop/IgorA', '-r300', '-dpng');

f                   = figure;
im                  = surf(Ra,Rb,RB);
im.LineStyle        = 'none';
ax                  = gca;
ax.Title.String     = 'FINITE: Reward Player B';
ax.XTick            = [0 .5 1];
ax.YTick            = [0 .5 1];
ax.ZTick            = [0 .5 1];
ax.FontSize         = 16;
ax.XLabel.String    = {'Score B'; '[Acc * Conf]'};
ax.YLabel.String    = {'Score A'; '[Acc * Conf]'};
ax.ZLabel.String    = 'Reward [a.u.]';
colormap turbo
print(f, '/Users/fschneider/Desktop/IgorB', '-r300', '-dpng');

%%% ZERO SUM %%%
for iA = 1:length(Ra)
    for iB = 1:length(Rb)
        if Ra(iA) >= Rb(iB)
            RA(iA,iB) 	= Ra(iA);
            RB(iA,iB) 	= -Ra(iA);
        else
            RA(iA,iB) 	= -Rb(iB);
            RB(iA,iB) 	= Rb(iB);
        end
    end
end

% Plot
f                   = figure;
im                  = surf(Ra,Rb,RA);
im.LineStyle        = 'none';
ax                  = gca;
ax.Title.String     = 'ZERO SUM: Reward Player A';
ax.XTick            = [0 .5 1];
ax.YTick            = [0 .5 1];
ax.ZTick            = [-1 -.5 0 .5 1];
ax.ZTickLabel       = {'-100','-50','0','50','100'};
ax.FontSize         = 16;
ax.XLabel.String    = {'Score B'; '[Acc * Conf]'};
ax.YLabel.String    = {'Score A'; '[Acc * Conf]'};
ax.ZLabel.String    = 'Reward [%]';
colormap turbo
print(f, '/Users/fschneider/Desktop/ZeroSumA', '-r300', '-dpng');

f                   = figure;
im                  = surf(Ra,Rb,RB);
im.LineStyle        = 'none';
ax                  = gca;
ax.Title.String     = 'ZERO SUM: Reward Player B';
ax.XTick            = [0 .5 1];
ax.YTick            = [0 .5 1];
ax.ZTick            = [-1 -.5 0 .5 1];
ax.ZTickLabel       = {'-100','-50','0','50','100'};
ax.FontSize         = 16;
ax.XLabel.String    = {'Score B'; '[Acc * Conf]'};
ax.YLabel.String    = {'Score A'; '[Acc * Conf]'};
ax.ZLabel.String    = 'Reward [%]';
colormap turbo
print(f, '/Users/fschneider/Desktop/ZeroSumB', '-r300', '-dpng');


%% CO-ACTION

for iA = 1:length(Ra)
    for iB = 1:length(Rb)
        Rsum(iA,iB) 	= sum([Ra(iA),Rb(iB)]) * R;
        Rmin(iA,iB) 	= min(Ra(iA),Rb(iB)) * R;
    end
end

% Plot
f                   = figure;
im                  = surf(Ra,Rb,Rmin);
im.LineStyle        = 'none';
ax                  = gca;
ax.Title.String     = 'COACTION: Min(RewA,RewB)';
ax.XTick            = [0 .5 1];
ax.YTick            = [0 .5 1];
ax.ZTick            = [0 .5 1];
ax.ZTickLabel       = {'0','50','100'};
ax.FontSize         = 16;
ax.XLabel.String    = {'Score B'; '[Acc * Conf]'};
ax.YLabel.String    = {'Score A'; '[Acc * Conf]'};
ax.ZLabel.String    = 'Reward [%]';
colormap turbo
print(f, '/Users/fschneider/Desktop/Coaction_min', '-r300', '-dpng');

f                   = figure;
im                  = surf(Ra,Rb,Rsum);
im.LineStyle        = 'none';
ax                  = gca;
ax.Title.String     = 'COACTION: Sum(RewA, RewB)';
ax.XTick            = [0 .5 1];
ax.YTick            = [0 .5 1];
ax.ZTick            = [0 .5 1 1.5 2];
ax.ZTickLabel       = {'0','50','100','150','200'};
ax.FontSize         = 16;
ax.XLabel.String    = {'Score B'; '[Acc * Conf]'};
ax.YLabel.String    = {'Score A'; '[Acc * Conf]'};
ax.ZLabel.String    = 'Reward [%]';
colormap turbo
print(f, '/Users/fschneider/Desktop/Coaction_sum', '-r300', '-dpng');

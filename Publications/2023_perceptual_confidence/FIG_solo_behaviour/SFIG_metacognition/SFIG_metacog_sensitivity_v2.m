close all

addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab

source_dir = '/Users/fschneider/ownCloud/var_plot/';
load([source_dir '/solo_performance.mat'])
load([source_dir '/hh_dyad_performance.mat'])
load([source_dir '/solo_correlation.mat'])
load([source_dir '/hh_dyad_correlation.mat'])

dest_dir            = '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG_solo_behaviour/SFIG_metacognition/raw/';
col                 = cool(7);
lw                  = 1;
lb_fs               = 8;

[x,y,auc]           = metacog(solo_perf, true, dest_dir);
[x_dy,y_dy,auc_dy]  = metacog(dyad_perf, false, dest_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check coherence effects %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = []; sub = []; c = [];
for i = 1:size(auc,1)
    data                = [data auc(i,:)];
    sub                 = [sub zeros(1,size(auc,2))+i];
    c                   = [c 1:size(auc,2)];
end
[P,T,STATS,TERMS]   = anovan(data,{sub,c},'varnames',{'Subject','Coherence'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Population AUC values %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f               	= figure('units','centimeters','position',[0 0 7.5 7.5]); hold on

vl                  = violinplot(auc);
vl                  = improveViolin(vl,col);
snr                 = solo_perf{1}.carr;
ax                 	= gca;
ax.XTick            = 1:length(vl);
ax.XTickLabel       = round(snr.*100);
ax.XLabel.String  	= 'Coherence';
ax.YLabel.String  	= 'AUC';
ax.YLim           	= [.35 1];
ax.FontSize         = lb_fs;
ax.Box              = 'off';

line([0 8],[.5 .5], 'color',[0 0 0], 'LineStyle', ':')
print(f, [dest_dir '/metacog_population'], '-r500', '-dpng');
print(f, [dest_dir '/metacog_population'], '-r500', '-dsvg');

% Test difference from 0.5
for iCoh = 1:size(auc,2)
    [Pval(iCoh),Hval(iCoh),STATS(iCoh)] = signrank(auc(:,iCoh),.5);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparison solo-dyadic %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f               	= figure('units','centimeters','position',[0 0 7.5 7.5]); hold on
auc_dy(auc_dy == 0) = nan;
vl                  = violinplot(auc - auc_dy);
vl                  = improveViolin(vl,col);

ax                 	= gca;
ax.XTick            = 1:length(vl);
ax.XTickLabel       = round(snr.*100);
ax.XLabel.String  	= 'Coherence';
ax.YLabel.String  	= 'AUC difference [Solo - Dyadic]';
ax.FontSize         = lb_fs;
ax.Box              = 'off';

line([0 8],[0 0], 'color',[0 0 0], 'LineStyle', ':')
print(f, [dest_dir '/metacog_difference'], '-r500', '-dpng');
print(f, [dest_dir '/metacog_difference'], '-r500', '-dsvg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example ROC %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f               	= figure('units','centimeters','position',[0 0 7.5 7.5]); hold on

ln                	= line([0 1],[0 1]);
ln.Color           	= [0 0 0];
ln.LineWidth       	= lw;
ln.LineStyle       	= ':';

for iCoh = 1:length(snr)
    plot(x{iSubj,iCoh},y{iSubj,iCoh},'Color',col(iCoh,:),'LineWidth',lw);
end

ax                 	= gca;
ax.YLabel.String  	= 'p(high eccentricity)';
ax.XLabel.String  	= 'p(low eccentricity)';
ax.XLim           	= [0 1];
ax.YLim           	= [0 1];
ax.FontSize         = lb_fs;
ax.Box              = 'off';

print(f, [dest_dir '/metacog_example_roc'], '-r500', '-dpng');
print(f, [dest_dir '/metacog_example_roc'], '-r500', '-dsvg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y,auc] = metacog(in_perf, plot_flag, dest_dir)

n_min_samples       = 5;
snr                 = in_perf{end}.carr;

for iSubj = 1:length(in_perf)
    
    if isempty(in_perf{iSubj})
        continue
    end
    
    clear ecc acc
    ecc             = in_perf{iSubj}.ecc_state; % eccentricity
    acc             = in_perf{iSubj}.acc_state; % accuracy
    
    for iCoh = 1:length(snr)
        acc_idx     = acc{iCoh} > median(acc{iCoh}); % Index: Above median accuracy
        ecc_idx     = ecc{iCoh} > median(ecc{iCoh}); % Index: Above median eccentricity
        
        dat1        = acc{iCoh}(ecc_idx); % accuracy(high eccentricity)
        dat2        = acc{iCoh}(~ecc_idx); % accuracy(low eccentricity)
        
        % Plot distributions
        if plot_flag && iSubj == 10
            if iCoh == 1
                f = figure('units','centimeters','position',[0 0 5 15]); hold on
            end
            plot_distributions(dat1,dat2,iCoh)
        end
        
        [x{iSubj,iCoh},y{iSubj,iCoh},auc(iSubj,iCoh)] = getAUROC(dat2, dat1);
    end
end

if plot_flag
    print(f, [dest_dir '/metacog_example_hist'], '-r500', '-dsvg', '-painters');
    print(f, [dest_dir '/metacog_example_hist'], '-r500', '-dpng', '-painters');
end
end

function [x,y,auc] = getAUROC(in1, in2)

if size(in1,1) == 1
    in1         = in1';
    in2         = in2';
end

lab          	= [zeros(length(in1),1); ones(length(in2),1)];
[x,y,~,auc]     = perfcurve(lab,[in1; in2],1);

end

function plot_distributions(dat1,dat2,iCoh)
ax                    	= subplot(7,1,iCoh); hold on
edges                   = 0:.025:1;
col1                    = cool(7);
col2                    = [.3 .3 .3];
alp                     = .5;
lb_fs                   = 8;

h1                     	= histogram(dat1,edges); hold on
h1.FaceColor          	= col1(iCoh,:);
h1.EdgeColor          	= 'none';
h1.FaceAlpha          	= alp;

h2                    	= histogram(dat2,edges);
h2.FaceColor           	= col2;
h2.EdgeColor          	= 'none';
h2.FaceAlpha          	= alp;

ax                     	= gca;
ax.XLim               	= [0 1];
ax.YLabel.String      	= '# States';
ax.XLabel.String      	= 'Accuracy [%]';
ax.XTick               	= 0:.25:1;
ax.FontSize           	= lb_fs;
ax.Box                	= 'off';

if iCoh < 7
    ax.XAxis.Visible        = 'off';
else
    legend('High','Low','Location','northwest')
end
end

function vl = improveViolin(vl,col_map)
for iV=1:length(vl)
    vl(iV).BoxWidth                     = .03;
    vl(iV).ViolinColor{1}               = col_map(iV,:);
    vl(iV).ViolinAlpha{1}               = .8;
    vl(iV).ScatterPlot.MarkerFaceColor  = [.25 .25 .25];
    vl(iV).ScatterPlot.MarkerEdgeColor  = 'none';
    vl(iV).ScatterPlot.SizeData         = 20;
end
end
close all

addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab

source_dir = '/Users/fschneider/ownCloud/var_plot/';
load([source_dir '/solo_performance.mat'])
load([source_dir '/hh_dyad_performance.mat'])

dest_dir            = '/Users/fschneider/Desktop/';
col                 = cool(7);
lw                  = 1;
lb_fs               = 8;

[x,y,auc]           = metacog(solo_perf);
[x_dy,y_dy,auc_dy]  = metacog(dyad_perf);

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

ax                 	= gca;
ax.XTick            = 1:length(snr);
ax.XTickLabel       = round(snr.*100);
ax.XLabel.String  	= 'Coherence';
ax.YLabel.String  	= 'AUC';
ax.YLim           	= [.35 1];
ax.FontSize         = lb_fs;
ax.Box              = 'off';
title('All states')

line([0 8],[.5 .5], 'color',[0 0 0], 'LineStyle', ':')
print(f, [dest_dir '/metacog_population'], '-r500', '-dpng');
print(f, [dest_dir '/metacog_population'], '-r500', '-dsvg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparison solo-dyadic %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f               	= figure('units','centimeters','position',[0 0 7.5 7.5]); hold on
auc_dy(auc_dy == 0) = nan;
vl                  = violinplot(auc - auc_dy);
vl                  = improveViolin(vl,col);

ax                 	= gca;
ax.XTick            = 1:length(snr);
ax.XTickLabel       = round(snr.*100);
ax.XLabel.String  	= 'Coherence';
ax.YLabel.String  	= 'AUC difference [Solo - Dyadic]';
ax.FontSize         = lb_fs;
ax.Box              = 'off';

line([0 8],[0 0], 'color',[0 0 0], 'LineStyle', ':')
print(f, [dest_dir '/metacog_difference'], '-r500', '-dpng');
print(f, [dest_dir '/metacog_difference'], '-r500', '-dsvg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example eccentricity distribution %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear ecc acc coh outc
iSubj         	= 10;
ecc             = solo_perf{iSubj}.ecc_state; % Target eccentricity
acc             = solo_perf{iSubj}.acc_state; % Target accuracy

f               = figure('units','centimeters','position',[0 0 5 15]); hold on
edges           = 0.5:.01:1;
col1            = [.3 .3 .3];
alp             = .5;

for iCoh = 1:length(snr)
    clear dat1 dat2 cIdx acc_idx
    
    ax          = subplot(7,1,iCoh);
    acc_idx     = acc{iCoh} > median(acc{iCoh}); % Index: Above median accuracy
    ecc_idx     = ecc{iCoh} > median(ecc{iCoh}); % Index: Above median eccentricity
    
    % Accuracy filtered for outcome and high eccentricity
%     dat1        = acc{iCoh}(ecc_idx); % p(accurate & high ecc)
%     dat2        = acc{iCoh}(~ecc_idx); % p(accurate & low ecc)    
    dat1        = acc{iCoh}(acc_idx & ecc_idx); % p(accurate & high ecc)
    dat2        = acc{iCoh}(acc_idx & ~ecc_idx); % p(accurate & low ecc)
    
    h1                          = histogram(dat1,edges); hold on
    h1.FaceColor                = col(iCoh,:);
    h1.EdgeColor                = 'none';
    h1.FaceAlpha                = alp;
    
    h2                          = histogram(dat2,edges);
    h2.FaceColor                = col1;
    h2.EdgeColor                = 'none';
    h2.FaceAlpha                = alp;
    
    ax                          = gca;
    ax.XLim                     = [.5 1];
    ax.YLabel.String            = '# States';
    ax.XLabel.String            = 'Accuracy [%]';
    
    ax.XTick                    = 0:.25:1;
    ax.FontSize                 = lb_fs;
    ax.Box                      = 'off';
    
    if iCoh < 7
        ax.XAxis.Visible        = 'off';
    end
    
    if iCoh == 7
        title('Filtered by eccentricity')
        legend('HighEcc','LowEcc','Location','northwest')
    end
end

print(f, [dest_dir '/metacog_example_hist'], '-r500', '-dsvg', '-painters');
print(f, [dest_dir '/metacog_example_hist'], '-r500', '-dpng', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example AUROC %%%
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
ax.YLabel.String  	= 'p(high accuracy & high eccentricity)'; %%% CHECK!
ax.XLabel.String  	= 'p(high accuracy & low eccentricity)';
ax.XLim           	= [0 1];
ax.YLim           	= [0 1];
ax.FontSize         = lb_fs;
ax.Box              = 'off';

print(f, [dest_dir '/metacog_example_roc'], '-r500', '-dpng');
print(f, [dest_dir '/metacog_example_roc'], '-r500', '-dsvg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y,auc] = metacog(in_perf)

n_min_samples       = 5;

for iSubj = 1:length(in_perf)
    
    if isempty(in_perf{iSubj})
        continue
    end
    
    clear ecc acc
    ecc             = in_perf{iSubj}.ecc_state; % Target eccentricity
    acc             = in_perf{iSubj}.acc_state; % Target accuracy
    
    for iCoh = 1:length(ecc)
        acc_idx     = acc{iCoh} > median(acc{iCoh}); % Index: Above median accuracy
        ecc_idx     = ecc{iCoh} > median(ecc{iCoh}); % Index: Above median eccentricity
       
        % Accuracy filtered for outcome and high eccentricity
        dat1        = acc{iCoh}(acc_idx & ecc_idx); % p(accurate & high ecc)
        dat2        = acc{iCoh}(acc_idx & ~ecc_idx); % p(accurate & low ecc)
%         dat1        = acc{iCoh}(ecc_idx); % p(accurate & high ecc)
%         dat2        = acc{iCoh}(~ecc_idx); % p(accurate & low ecc)
        
        % Filter for number of samples
        if  length(dat1) >= n_min_samples && length(dat2) >= n_min_samples
            [x{iSubj,iCoh},y{iSubj,iCoh},auc(iSubj,iCoh)] = getAUROC(dat2, dat1);
        else
            x{iSubj,iCoh}   = nan;
            y{iSubj,iCoh}   = nan;
            auc(iSubj,iCoh) = nan;
        end
    end
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

close all

addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab

source_dir = '/Users/fschneider/ownCloud/var_plot/';
load([source_dir '/solo_performance.mat'])
load([source_dir '/hh_dyad_performance.mat'])
load([source_dir '/solo_correlation.mat'])
load([source_dir '/hh_dyad_correlation.mat'])

dest_dir            = '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG_metacogntion/raw/';
col                 = cool(7);
lw                  = 1;
lb_fs               = 8;
option_flag         = 1;

[x,y,auc]           = metacog(solo_perf,solo_cr,option_flag);
[x_dy,y_dy,auc_dy]  = metacog(dyad_perf,solo_cr,option_flag);

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

line([0 8],[.5 .5], 'color',[0 0 0], 'LineStyle', ':')
print(f, [dest_dir '/metacog_population_'  num2str(option_flag)], '-r500', '-dpng');
print(f, [dest_dir '/metacog_population_'  num2str(option_flag)], '-r500', '-dsvg');

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
print(f, [dest_dir '/metacog_difference_'  num2str(option_flag)], '-r500', '-dpng');
print(f, [dest_dir '/metacog_difference_'  num2str(option_flag)], '-r500', '-dsvg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example eccentricity distribution %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear ecc acc coh outc
iSubj         	= 10;
ecc             = solo_perf{iSubj}.trg_all_ecc; % Target eccentricity
acc             = solo_perf{iSubj}.trg_all_acc; % Target accuracy
coh             = solo_perf{iSubj}.trg_all_coh; % Target coherence
outc            = solo_perf{iSubj}.trg_all_outc; % Target outcome
ts              = solo_perf{iSubj}.trg_ts_state; % Target timestamp

f               = figure('units','centimeters','position',[0 0 5 15]); hold on
edges           = 0:.05:1;
col1            = [.3 .3 .3];
alp             = .5;

for iCoh = 1:length(snr)
    clear dat1 dat2 cIdx acc_idx
    
    ax          = subplot(7,1,iCoh);
    cIdx        = coh == snr(iCoh); % Coherence index
    acc_idx     = acc > median(acc(cIdx)); % Index: Above median accuracy
    ecc_idx     = ecc > median(ecc(cIdx)); % Index: Above median eccentricity
    ts_idx      = ts > 500; % Index: Target presentation >1000ms after direction change
    
    if option_flag == 1
        % Option 1: eccentricity filtered for hits and accuracy
        dat1        = ecc(cIdx & ts_idx & outc & acc_idx); % p(correct & high accuracy)
        dat2        = ecc(cIdx & ts_idx & outc & ~acc_idx); % p(correct & low accuracy)
    elseif option_flag == 2
        % Option 2: accuracy filtered for outcome and high eccentricity
        dat1        = acc(cIdx & ts_idx & outc & ecc_idx); % p(correct & high ecc)
        dat2        = acc(cIdx & ts_idx & ~outc & ecc_idx); % p(incorrect & high ecc)
    end
    
    h1                          = histogram(dat1,edges); hold on
    h1.FaceColor                = col1;
    h1.EdgeColor                = 'none';
    h1.FaceAlpha                = alp;
    
    h2                          = histogram(dat2,edges);
    h2.FaceColor                = col(iCoh,:);
    h2.EdgeColor                = 'none';
    h2.FaceAlpha                = alp;
    
    ax                          = gca;
    ax.XLim                     = [0 1];
    ax.YLabel.String            = '# Targets';
    
    if option_flag == 1
        ax.XLabel.String      	= 'Eccentricity [%]';
    elseif option_flag == 2
        ax.XLabel.String      	= 'Accuracy [%]';
    end
    
    ax.XTick                    = 0:.25:1;
    ax.FontSize                 = lb_fs;
    ax.Box                      = 'off';
    
    if iCoh < 7
        ax.XAxis.Visible        = 'off';
    end
    
    if iCoh == 1
        if option_flag == 1
            title('Filtered by eccentricity')
        elseif option_flag == 2
            title('Filtered by accuracy')
        end
        legend('HI','MI','Location','northeast')
    end
end

print(f, [dest_dir '/metacog_example_hist_' num2str(option_flag)], '-r500', '-dsvg', '-painters');
print(f, [dest_dir '/metacog_example_hist_' num2str(option_flag)], '-r500', '-dpng', '-painters');

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
if option_flag == 1
    ax.YLabel.String  	= 'p(high accuracy & correct)';
    ax.XLabel.String  	= 'p(low accuracy & correct)';
elseif option_flag == 2
    ax.YLabel.String  	= 'p(high eccentricity & incorrect)';
    ax.XLabel.String  	= 'p(high eccentricity & correct)';
end

ax.XLim           	= [0 1];
ax.YLim           	= [0 1];
ax.FontSize         = lb_fs;
ax.Box              = 'off';

print(f, [dest_dir '/metacog_example_roc_' num2str(option_flag)], '-r500', '-dpng');
print(f, [dest_dir '/metacog_example_roc_' num2str(option_flag)], '-r500', '-dsvg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [x,y,auc] = metacog_state(in_perf,in_cr,option_flag)
% 
% n_min_samples       = 5;
% snr                 = in_perf{end}.carr;
% 
% for iSubj = 1:length(in_perf)
%     
%     if isempty(in_perf{iSubj})
%         continue
%     end
%     
%     clear ecc acc coh outc ts
%     ecc             = in_perf{iSubj}.trg_all_ecc; % Target eccentricity
%     acc             = in_perf{iSubj}.trg_all_acc; % Target accuracy
%     coh             = in_perf{iSubj}.trg_all_coh; % Target coherence
%     
%     for iCoh = 1:length(snr)
%         cIdx        = coh == snr(iCoh); % Coherence index
%         acc_idx     = acc > median(acc(cIdx)); % Index: Above median accuracy
%         ecc_idx     = ecc > median(ecc(cIdx)); % Index: Above median eccentricity
%         ts_idx      = ts > mean(in_cr{iSubj}.lag); % Index: Target presentations after average response lag
%         
%         if option_flag == 1
%             % Option 1: eccentricity filtered for hits and accuracy
%             dat1        = ecc(cIdx & ts_idx & outc & acc_idx); % p(correct & high accuracy)
%             dat2        = ecc(cIdx & ts_idx & outc & ~acc_idx); % p(correct & low accuracy)
%         elseif option_flag == 2
%             % Option 2: accuracy filtered for outcome and high eccentricity
%             dat1        = acc(cIdx & ts_idx & outc & ecc_idx); % p(correct & high ecc)
%             dat2        = acc(cIdx & ts_idx & ~outc & ecc_idx); % p(incorrect & high ecc)
%         end
%         
%         % Filter for number of samples
%         if  length(dat1) >= n_min_samples && length(dat2) >= n_min_samples
%             [x{iSubj,iCoh},y{iSubj,iCoh},auc(iSubj,iCoh)] = getAUROC(dat2, dat1);
%         else
%             x{iSubj,iCoh}   = nan;
%             y{iSubj,iCoh}   = nan;
%             auc(iSubj,iCoh) = nan;
%         end
%     end
% end
% end


function [x,y,auc] = metacog(in_perf,in_cr,option_flag)

n_min_samples       = 5;
snr                 = in_perf{end}.carr;

for iSubj = 1:length(in_perf)
    
    if isempty(in_perf{iSubj})
        continue
    end
    
    clear ecc acc coh outc ts
    ecc             = in_perf{iSubj}.trg_all_ecc; % Target eccentricity
    acc             = in_perf{iSubj}.trg_all_acc; % Target accuracy
    coh             = in_perf{iSubj}.trg_all_coh; % Target coherence
    outc            = in_perf{iSubj}.trg_all_outc; % Target outcome
    ts              = in_perf{iSubj}.trg_ts_state; % Target timestamp
    
    for iCoh = 1:length(snr)
        cIdx        = coh == snr(iCoh); % Coherence index
        acc_idx     = acc > median(acc(cIdx)); % Index: Above median accuracy
        ecc_idx     = ecc > median(ecc(cIdx)); % Index: Above median eccentricity
        ts_idx      = ts > mean(in_cr{iSubj}.lag); % Index: Target presentations after average response lag
        
        if option_flag == 1
            % Option 1: eccentricity filtered for hits and accuracy
            dat1        = ecc(cIdx & ts_idx & outc & acc_idx); % p(correct & high accuracy)
            dat2        = ecc(cIdx & ts_idx & outc & ~acc_idx); % p(correct & low accuracy)
        elseif option_flag == 2
            % Option 2: accuracy filtered for outcome and high eccentricity
            dat1        = acc(cIdx & ts_idx & outc & ecc_idx); % p(correct & high ecc)
            dat2        = acc(cIdx & ts_idx & ~outc & ecc_idx); % p(incorrect & high ecc)
        end
        
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
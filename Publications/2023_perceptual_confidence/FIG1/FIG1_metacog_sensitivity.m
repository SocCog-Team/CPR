addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/solo_performance.mat')

dest_dir            = '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG1/raw';
snr                 = solo_perf{end}.carr;
col                 = cool(length(snr));
lw                  = 1;
lb_fs               = 8;

% Data processing
for iSubj = 1:length(solo_perf)
    clear ecc acc coh outc
    ecc             = solo_perf{iSubj}.trg_all_ecc; % Target eccentricity
    acc             = solo_perf{iSubj}.trg_all_acc; % Target accuracy
    coh             = solo_perf{iSubj}.trg_all_coh; % Target coherence
    outc            = solo_perf{iSubj}.trg_all_outc; % Target outcome
    
    for iCoh = 1:length(snr)
        cIdx        = coh == snr(iCoh); % Coherence index
        acc_idx     = acc > median(acc(cIdx)); % Index: Above median accuracy
        ecc_idx     = ecc > median(ecc(cIdx)); % Index: Above median eccentricity
        
%         % Option 1: filtered for hits and accuracy   
%         dat1        = ecc(cIdx & outc & acc_idx); % p(correct & high accuracy)
%         dat2        = ecc(cIdx & outc & ~acc_idx); % p(correct & low accuracy)
                
        % Option 1: filtered for hits and eccentricity
        dat1        = acc(cIdx & outc & ecc_idx); % p(correct & high ecc)
        dat2        = acc(cIdx & ~outc & ecc_idx); % p(incorrect & high ecc)

                
        % Option 2: filtered for high eccentricity and outcome   
%         dat1 = ecc(cIdx & ~ecc_idx & outc); % p(high conf & correct)
%         dat2 = ecc(cIdx & ~ecc_idx & ~outc); % p(high conf & incorrect)
                
%         % Option 3: filtered for high eccentricity and accuracy   
%         dat1 = ecc(cIdx & ecc_idx & acc_idx); % p(high conf & high accuracy)
%         dat2 = ecc(cIdx & ecc_idx & ~acc_idx); % p(high conf & low accuracy)
        
        [x{iSubj,iCoh},y{iSubj,iCoh},auc(iSubj,iCoh)] = getAUROC(dat2, dat1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example eccentricity distribution %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear ecc acc coh outc
iSubj         	= 10;
ecc             = solo_perf{iSubj}.trg_all_ecc;
acc             = solo_perf{iSubj}.trg_all_acc;
coh             = solo_perf{iSubj}.trg_all_coh;
outc            = solo_perf{iSubj}.trg_all_outc;

f               = figure('units','centimeters','position',[0 0 5 15]); hold on
edges           = 0:.05:1;
col1            = [.3 .3 .3];
alp             = .5;

for iCoh = 1:length(snr)
    clear dat1 dat2 cIdx acc_idx
    
    ax          = subplot(7,1,iCoh);
    cIdx        = coh == snr(iCoh);
    acc_idx     = acc > median(acc(cIdx));
    ecc_idx     = ecc > median(ecc(cIdx)); % Above median eccentricity
    
%     dat1        = ecc(cIdx & outc & acc_idx); % p(correct & high accuracy)
%     dat2        = ecc(cIdx & outc & ~acc_idx); % p(correct & low accuracy)

dat1        = acc(cIdx & outc & ecc_idx); % p(correct & high ecc)
dat2        = acc(cIdx & ~outc & ecc_idx); % p(incorrect & high ecc)

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
    ax.XLabel.String            = 'Accuracy [%]';
    ax.XTick                    = 0:.25:1;
    ax.FontSize                 = lb_fs;
    ax.Box                      = 'off';
    
    if iCoh < 7
        ax.XAxis.Visible        = 'off';
    end
    
    if iCoh == 1
        title('Filtered by accuracy')
        legend('HI','MI','Location','northeast')
    end
end

% lg                          = legend('High Accuracy', 'Low Accuracy', 'Location', 'eastoutside');
% box off
% axis square

print(f, [dest_dir '/metacog_example_hist'], '-r500', '-dsvg', '-painters');
print(f, [dest_dir '/metacog_example_hist'], '-r500', '-dpng', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example AUROC %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f               	= figure('units','centimeters','position',[0 0 7.5 7.5]); hold on
ex_subj             = 10;

ln                	= line([0 1],[0 1]);
ln.Color           	= [0 0 0];
ln.LineWidth       	= lw;
ln.LineStyle       	= ':';

for iCoh = 1:length(snr)
    plot(x{ex_subj,iCoh},y{ex_subj,iCoh},'Color',col(iCoh,:),'LineWidth',lw);
end

ax                 	= gca;
ax.YLabel.String  	= 'p(high eccentricity & incorrect)';
ax.XLabel.String  	= 'p(high eccentricity & correct)';
ax.XLim           	= [0 1];
ax.YLim           	= [0 1];
ax.FontSize         = lb_fs;
ax.Box              = 'off';

print(f, [dest_dir '/metacog_example_roc'], '-r500', '-dpng');
print(f, [dest_dir '/metacog_example_roc'], '-r500', '-dsvg');

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
ax.YLim           	= [.45 1];
ax.FontSize         = lb_fs;
ax.Box              = 'off';

line([0 8],[.5 .5], 'color',[0 0 0], 'LineStyle', ':')
print(f, [dest_dir '/metacog_population'], '-r500', '-dpng');
print(f, [dest_dir '/metacog_population'], '-r500', '-dsvg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
close all
clear all

addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab

% Adjust path
source_pth = '/Users/fschneider/Documents/GitHub/CPR/Publications/2024_perceptual_confidence/var_plot/';
load([ source_pth '/solo_correlation.mat'])
load([ source_pth '/solo_performance.mat'])

lw                  = 1.5;
lb_fs               = 8;
alp                 = .2;
frme_ms             = 1000/120;
nLag                = 150;
snr                 = unique(solo_cr{1}.coh);

%% Example subject
subj_id             = 17;
col                 = cool(length(snr));

f                	= figure('units','centimeters','position',[0 0 6 5]);
hold on

for iCoh = 1:length(snr)
    idx             = solo_cr{subj_id}.coh == snr(iCoh);
    dat             = solo_cr{subj_id}.sxc(idx,:);
    coh_id{iCoh}   	= num2str(round(snr(iCoh),2));
        
    pk           	= find(mean(dat(:,150:end)) == max(mean(dat(:,150:end))));
    ln            	= line([pk pk],[0 max(mean(dat(:,150:end)))], 'Color', col(iCoh, :),'LineStyle', ':', 'LineWidth', lw/1.5);
   
    % Boostrap confidence intervals
    nRep        	= 1000;
    [CI,~]        	= bootci(nRep,{@mean,dat(:,150:end)},'Alpha',0.05);
    
    % Prepare filled area
    vec           	= 1:length(CI);
    x_spacing     	= [vec fliplr(vec)];
    ci            	= [CI(1,:) fliplr(CI(2,:))];

    % Overlay confidence intervals
    fl           	= fill(x_spacing,ci,col(iCoh, :),'EdgeColor','none', 'FaceAlpha', alp);
    
    % Plot mean curve
    pl(iCoh)        = plot(mean(dat(:,150:end)), 'LineWidth', lw, 'Color', col(iCoh, :)); 
end

ax                  = gca;
ax.XLabel.String    = 'Lag [ms]';
ax.XLim             = [1 150];
ax.XTick            = [1 60 120];
ax.XTickLabel       = round(cellfun(@str2num, ax.XTickLabel) * frme_ms);
ax.YLabel.String    = 'XCorr coef [norm]';
ax.YLim             = [0 .15];
ax.FontSize      	= lb_fs;

% lg                  = legend(pl,coh_id);
% lg.FontSize         = lb_fs;
% lg.Box              = 'off';
% lg.Location         = 'northwest';

dest_dir            = '/Users/fschneider/Documents/GitHub/CPR/Publications/2024_perceptual_confidence/FIG_solo_behaviour/raw/';
print(f, [dest_dir '/subj_correlation'], '-r500', '-dpng');
print(f, [dest_dir '/subj_correlation'], '-r500', '-dsvg');

%% Peak quantification

for iSubj = 1:length(solo_cr)
    for iCoh = 1:length(snr)
        idx                 = solo_cr{iSubj}.coh == snr(iCoh);
        dat                 = solo_cr{iSubj}.sxc(idx,:);

        for iBlock = 1:size(dat,1)
            % Find peak lag
            pk_pos{iSubj,iCoh}(iBlock) = find(max(dat(iBlock,150:end)) == dat(iBlock,150:end));
        end
    end
end

avg_pk_pos = cellfun(@median,pk_pos);
std_pk_pos = cellfun(@std,pk_pos);

f                   = figure('units','centimeters','position',[0 0 5 6]); hold on
ax = subplot(2,1,1);
vl                  = violinplot( (avg_pk_pos * frme_ms) ./1e3);
vl                  = improveViolin(vl,col);
ax                  = gca;
ax.YLabel.String    = 'Median lag [s]';
ax.XTick            = 1:7;
ax.XTickLabel       = round(snr,2)*100;
ax.XLabel.String    = 'Coherence [%]';
ax.FontSize      	= lb_fs;
ax.Box              = 'off';
ax.XAxis.Visible    = 'off';
axis tight

ax                  = subplot(2,1,2);
vl                  = violinplot(((std_pk_pos) * frme_ms) ./ 1e3);
vl                  = improveViolin(vl,col);
ax                  = gca;
ax.YLabel.String    = 'Std lag [s]';
ax.XTick            = 1:7;
ax.XTickLabel       = round(snr,2)*100;
ax.XLabel.String    = 'Coherence [%]';
ax.FontSize      	= lb_fs;
ax.Box              = 'off';
axis tight

print(f, [dest_dir '/pop_lag'], '-r500', '-dpng');
print(f, [dest_dir '/pop_lag'], '-r500', '-dsvg');

% bx                  = boxplot((avg_pk_pos - nLag) * frme_ms, 'Color', [0 0 0]);
% set(bx,'MarkerEdgeColor','k')
% set(bx,'LineWidth', 2);
% h = findobj(ax,'Tag','Box');
% for j=1:length(h)
%     patch(get(h(j),'XData'),get(h(j),'YData'),col(j,:),'FaceAlpha',.5);
% end 

%% Peak-Average Ratio
 
% for iSub = 1:size(solo_cr,2)
%     
%     sxc = solo_cr{iSub}.sxc;
%     coh = solo_cr{iSub}.coh;
%     for iCoh = 1:size(snr,2)
%         
%         msxc = mean(sxc(coh == snr(iCoh),:));
%         par(iSub,iCoh)   	= max(msxc) / mean(msxc);
%     end
% end
% 
% f                           = figure('units','centimeters','position',[0 0 6 5]);
% vl                          = violinplot(par);
% vl                        	= improveViolin(vl,col);
% ax                          = gca;
% ax.YLabel.String            = 'Response reliability';
% ax.XTick                    = 1:7;
% ax.XLim                     = [.5 7.5];
% ax.XTickLabel               = round(snr,2);
% ax.XLabel.String            = 'Coherence [%]';
% ax.FontSize                 = lb_fs;
% 
% print(f, [dest_dir '/pop_response_reliability'], '-r500', '-dpng');
% print(f, [dest_dir '/pop_response_reliability'], '-r500', '-dsvg');


%% FUNCTIONS %%%

function vl = improveViolin(vl,col_map)
for iV=1:length(vl)
    vl(iV).BoxWidth                     = .025;
    vl(iV).ViolinColor{1}               = col_map(iV,:);
    vl(iV).ViolinAlpha{1}               = .8;
    vl(iV).ScatterPlot.MarkerFaceColor  = [.25 .25 .25];
    vl(iV).ScatterPlot.MarkerEdgeColor  = 'none';
    vl(iV).ScatterPlot.SizeData         = 10;
end
end
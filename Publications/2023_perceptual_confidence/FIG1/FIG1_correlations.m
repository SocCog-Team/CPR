%%% Run FIG2 script first to get solo correlation data [cell, solo_cr] %%%

lw                  = 3;
lb_fs               = 16;
alp                 = .2;
frme_ms             = 1000/120;
nLag                = 150;
subj_id             = 17;
snr                 = unique(solo_cr{1}.coh);
col                 = cool(length(snr));

f = figure;hold on
ln                  = line([150 150],[0 .15], 'LineWidth', lw, 'LineStyle', ':', 'Color', [0 0 0]);

for iCoh = 1:length(snr)
    idx             = solo_cr{subj_id}.coh == snr(iCoh);
    dat             = solo_cr{subj_id}.sxc(idx,:);
    coh_id{iCoh}   	= num2str(round(snr(iCoh),2));
    
    % Boostrap confidence intervals
    nRep                        = 1000;
    [CI,~]                      = bootci(nRep,{@mean,dat},'Alpha',0.05);
    
    % Prepare filled area
    vec                         = 1:length(CI);
    x_spacing                   = [vec fliplr(vec)];
    ci                          = [CI(1,:) fliplr(CI(2,:))];

    % Overlay confidence intervals
    fl                          = fill(x_spacing,ci,col(iCoh, :),'EdgeColor','none', 'FaceAlpha', alp);
    
    % Plot mean curve
    pl(iCoh)                 	= plot(mean(dat), 'LineWidth', lw, 'Color', col(iCoh, :)); 
end

ax                  = gca;
ax.XLabel.String    = 'Lag [ms]';
ax.XLim             = [0 301];
ax.XTick            = [0 75 150 225 300];
ax.XTickLabel       = round((cellfun(@str2num, ax.XTickLabel)-nLag) * frme_ms);
ax.YLabel.String    = 'XCorr coef [norm]';
ax.YLim             = [0 .15];
ax.FontSize      	= lb_fs;

lg                  = legend(pl,coh_id);
lg.FontSize         = lb_fs;
lg.Box              = 'off';
lg.Location         = 'northwest';

dest_dir            = '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG1/raw';
print(f, [dest_dir '/subj_correlation'], '-r500', '-dpng');
print(f, [dest_dir '/subj_correlation'], '-r500', '-dsvg');

%% Peak quantification

for iSubj = 1:length(solo_cr)
    for iCoh = 1:length(snr)
        idx                 = solo_cr{iSubj}.coh == snr(iCoh);
        dat                 = solo_cr{iSubj}.sxc(idx,:);
        
        for iBlock = 1:size(dat,1)
            % Find peak lag
            pk_pos{iSubj,iCoh}(iBlock) = find(max(dat(iBlock,150:end)) == dat(iBlock,:));
        end
    end
end

avg_pk_pos = cellfun(@median,pk_pos);
std_pk_pos = cellfun(@std,pk_pos);

f                   = figure;
bx                  = boxplot((avg_pk_pos - nLag) * frme_ms, 'Color', [0 0 0]);
ax                  = gca;
ax.YLabel.String    = 'Avg lag [ms]';
ax.XTick            = 1:7;
ax.XTickLabel       = round(snr,2);
ax.XLabel.String    = 'Coherence [%]';
ax.FontSize      	= lb_fs;
ax.Box              = 'off';
set(bx,'LineWidth', 2);
print(f, [dest_dir '/pop_lag_median'], '-r500', '-dpng');
print(f, [dest_dir '/pop_lag_median'], '-r500', '-dsvg');

f                   = figure;
bx                  = boxplot(std_pk_pos, 'Color', [0 0 0]);
ax                  = gca;
ax.YLabel.String    = 'Std lag [ms]';
ax.XTick            = 1:7;
ax.XTickLabel       = round(snr,2);
ax.XLabel.String    = 'Coherence [%]';
ax.FontSize      	= lb_fs;
ax.Box              = 'off';
set(bx,'LineWidth', 2);
print(f, [dest_dir '/pop_lag_std'], '-r500', '-dpng');
print(f, [dest_dir '/pop_lag_std'], '-r500', '-dsvg');

%% Peak-Average Ratio

f                           = figure;hold on
msxc                     	= cellfun(@mean,sxc,'UniformOutput',false);

for iSub = 1:size(sxc,1)
    for iCoh = 1:size(sxc,2)
        x                	= msxc{iSub,iCoh};
        par(iSub,iCoh)   	= max(x) / mean(x);
    end
end

ax                         = gca;
[ax,pl]                    = plotData(ax,par,false,lw,alp,exmpl_id,exmpl_col,avg_mult);
ax.XLim                    = [1 7];
ax.XTick                   = [1 4 7];
ax.XTickLabel              = round([snr(1) snr(4) snr(7)],2)*100;
ax.XLabel.String           = {'Coherence'; '[%]'};
ax.YLabel.String           = {'Peak-Average ratio'};
ax.FontSize                = lb_fs;
ax.XTickLabelRotation      = 0;

%% Response lag - RT comparison

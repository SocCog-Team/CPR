f                           = figure('units','normalized','position',[0 0 .4 .635]);
lb_fs                       = 20;
lw                          = 3;
cmap                        = [0 0 0; gray(256)];
cmap                        = flipud(cmap);
pth                         = '/Volumes/T7_Shield/CPR_psychophysics/';

for iSubj = 1:2
    
    if iSubj == 1
        tbl             	= load([pth '/Dyad22/summary/20220510_anm_CPRdyadic_block1_tbl.mat']);
    else
        tbl             	= load([pth '/Dyad22/summary/20220510_akn_CPRdyadic_block1_tbl.mat']);
    end
    
    c                       = 0;
    t                    	= tbl.t;
    
    for iState = 1:size(t,1)
        if t.trg_shown(iState) == false
            continue
        end
        
        for iTarget = 1:length(t.trg_ts{iState})
            c                   = c+1;
            trg_acc{iSubj}(c)  	= t.trg_acc{iState}(iTarget);
            trg_conf{iSubj}(c)	= t.trg_ecc{iState}(iTarget);
        end
    end
end

% Calculate reward matrix
acc                         = 0:.001:1;
conf                        = 0:.001:1;
rew                         = acc' .* conf;

% Determine arc width for each confidence level
for j = 1:length(conf)
    arc(j)                  = 180 - (180 * conf(j));
end

% Cap arc width at target width (2dva == 12.7587deg at chosen position)
trg_idx                  	= arc < 12.7587;
arc(trg_idx)                = 12.7587;

% For each confidence level, calculate minimum accuracy required to hit
% the target at given arc width - normalised values
hit_width_acc               = 1 - ((arc/2) / 180);
hit_width_acc(trg_idx)   	= 1 - (12.7587/2)/180; % arc width fixed

% Remove position that cannot yield reward from reward matrix
for iAcc = 1:length(acc)
    indx                    = conf < hit_width_acc(iAcc);
    rew(iAcc,indx)          = nan;
end

% Plot reward matrix
ax0                         = gca; hold on
im                          = imagesc(acc,conf,rew);
ax0.XLabel.String           = 'Accuracy';
ax0.YLabel.String           = 'Eccentricity';
ax0.FontSize                = lb_fs;
ax0.XLim                    = [0 1];
ax0.YLim                    = [0 1];
ax0.XLim                    = [0 1];
ax0.XTick                   = [0:.2:1];
ax0.YTick                   = [0:.2:1];
ax0.XTickLabelRotation      = 0;

colormap(cmap)
cb                          = colorbar;
cb.Label.String             = 'Score';
cb.Location                 = 'northoutside';

% Add targets
sc1                         = scatter(trg_acc{1}, trg_conf{1}, 'filled');
sc1.CData                   = [.65 0 0];
sc1.SizeData                = 40;
sc1.MarkerFaceAlpha         = .75;

sc2                         = scatter(trg_acc{2}, trg_conf{2}, 'filled');
sc2.CData                   = [0 .65 0];
sc2.SizeData                = 40;
sc2.MarkerFaceAlpha         = .75;

ax0.Position                = [.25 .19 .65 .65];

% Add histograms on both axes
nBin                        = 25;
vofs                        = .175;
hofs                        = .175;

axh                        = axes('Position', [ax0.Position(1)-vofs ax0.Position(2) ax0.Position(3)/10 ax0.Position(4)]); hold on
for iSubj = 1:2
    [h, edg]            	= histcounts(trg_conf{iSubj},nBin);
    cntr                 	= edg(1:end-1) + diff(edg) ./ 2;
    st                     	= stairs(-h,cntr);
    st.LineWidth          	= lw;
    if iSubj == 1
        st.Color          	= [.65 0 0];
    else
        st.Color          	= [0 .65 0];
    end
    axh.YLim              	= [0 1];
    axh.XAxis.Visible     	= 'off';
    axh.YAxis.Visible       = 'off';
end

axv                         = axes('Position', [ax0.Position(1) ax0.Position(2)-vofs ax0.Position(3) ax0.Position(4)/10]); hold on
for iSubj = 1:2
    [v, edg]              	= histcounts(trg_acc{iSubj},nBin);
    cntr                   	= edg(1:end-1) + diff(edg) ./ 2;
    st                     	= stairs(cntr,-v);
    st.LineWidth          	= lw;
    if iSubj == 1
        st.Color          	= [.65 0 0];
    else
        st.Color          	= [0 .65 0];
    end
    axv.XLim              	= [0 1];
    axv.XAxis.Visible     	= 'off';
    axv.YAxis.Visible    	= 'off';
end

dest_dir                    = '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG_task/raw/';

print(f, [dest_dir '/target_distribution_dyad'], '-r500', '-dpng');
print(f, [dest_dir '/target_distribution_dyad'], '-r500', '-dsvg');

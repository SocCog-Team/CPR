f                           = figure('units','normalized','position',[0 0 .4 .5]);
lb_fs                       = 16;
lw                          = 3;
cmap                        = [0 0 0; gray(256)];
c                           = 0;
pth                         = '/Volumes/T7_Shield/CPR_psychophysics/';
tbl                         = load([pth 'AnM/summary/20220629_anm_CPRsolo_block2_tbl.mat']);
exmpl_col                   = [.5 0 .5 .5];
t                           = tbl.t;

% Extract experimental data
clear trg_acc trg_conf trg_coh trg_hit
for iState = 1:size(t,1)
    if t.trg_shown(iState) == false
        continue
    end
    
    for iTarget = 1:length(t.trg_ts{iState})
        c                   = c+1;
        trg_acc(c)          = t.trg_acc{iState}(iTarget);
        trg_conf(c)         = t.trg_ecc{iState}(iTarget);
        trg_coh(c)          = t.rdp_coh(iState);
        trg_hit(c)          = t.trg_hit{iState}(iTarget);
        
        if trg_acc(c) < .5 && trg_hit(c) == 1
            disp([num2str(iState) ' ' num2str(iTarget) ' ' num2str(c) ])
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
aidx                        = arc < 12.7587;
arc(aidx)                   = 12.7587;

% For each confidence level, calculate minimum accuracy required to hit
% the target at given arc width - normalised values
hit_width_acc               = 1 - ((arc/2) / 180);
hit_width_acc(aidx)         = 1 - (12.7587/2)/180; % arc width fixed

% Remove position that cannot yield reward from reward matrix
for iAcc = 1:length(acc)
    indx                    = conf < hit_width_acc(iAcc);
    rew(iAcc,indx)          = nan;
end

% Plot reward matrix
ax0 = gca;
hold on
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
cb.Label.String             = '% Reward';
cb.Location                 = 'eastoutside';

% Adjust axes position
ax0.Position(1:2)           = [.15 .17];
hofs                        = .135;
vofs                        = .16;
snr                         = unique(t.rdp_coh);

% cmap_coh                    = jet(size(snr,1));
cmap_coh                    = cool(size(snr,1));

% Add targets
for iCoh = 1:length(snr)
    cidx                    = trg_coh == snr(iCoh);
    sc                      = scatter(trg_acc(cidx), trg_conf(cidx), 'filled');
    sc.CData            	= cmap_coh(iCoh,:);
    sc.SizeData             = 50;
    sc.MarkerFaceAlpha      = .85;
end

% Add histograms on both axes
nBin                        = 25;
ax0h                        = axes('Position', [ax0.Position(1)-hofs ax0.Position(2) ax0.Position(3)/15 ax0.Position(4)]); hold on
[h, edg]                    = histcounts(trg_conf,nBin);
cntr                        = edg(1:end-1) + diff(edg) ./ 2;
st                          = stairs(-h,cntr);
st.LineWidth                = lw/1.5;
st.Color                    = [0 0 0];
ax0h.YLim                   = [0 1];
ax0h.XAxis.Visible          = 'off';
ax0h.YAxis.Visible          = 'off';

ax0v                        = axes('Position', [ax0.Position(1) ax0.Position(2)-vofs ax0.Position(3) ax0.Position(4)/15]); hold on
[v, edg]                    = histcounts(trg_acc,nBin);
cntr                        = edg(1:end-1) + diff(edg) ./ 2;
st                          = stairs(cntr,-v);
st.LineWidth                = lw/1.5;
st.Color                    = [0 0 0];
ax0v.XLim                   = [0 1];
ax0v.XAxis.Visible          = 'off';
ax0v.YAxis.Visible          = 'off';

dest_dir                    = '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG1/raw';
print(f, [dest_dir '/target_distribution'], '-r500', '-dpng');
print(f, [dest_dir '/target_distribution'], '-r500', '-dsvg');

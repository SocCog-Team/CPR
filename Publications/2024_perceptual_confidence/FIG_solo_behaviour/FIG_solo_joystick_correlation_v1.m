%%% Crosscorrelation between vector length and joystick displacement

addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
 
nLag        = 150;
smooth_win  = 20;
n_ofs       = 3;
nshift      = 30;

clear avg_xc_acc_vec avg_xc_ecc_vec avg_xc_ecc_acc avg_xc_control
    
load('/Users/fschneider/Desktop/resultant_vec/subj_lst.mat')

for iSubj = 1:length(sbj_lst)
    
    clear frme_vec js_ecc js_dir xc_acc_vec xc_ecc_vec xc_ecc_acc state_coh xc_control
    load(['/Users/fschneider/Desktop/resultant_vec/vector_subj_' sbj_lst{iSubj} '.mat'])
    load(['/Users/fschneider/Desktop/resultant_vec/joystick_subj_' sbj_lst{iSubj} '.mat'])
    cnt = 0;
    
    for iState = 1:size(frme_vec,2)
        clear v_smooth v_smooth_detrend ecc_detrend js_dir_error joystick_acc acc_detrend
        
        if length(js_ecc{iState}) == length(frme_vec{iState}.resultant_length)+1
            cnt = cnt+1;
            v_smooth = smoothdata(frme_vec{iState}.resultant_length,'gaussian',smooth_win);
            v_smooth_detrend = [0 v_smooth - nanmean(v_smooth)];
            
            ecc_detrend = js_ecc{iState} - nanmean(js_ecc{iState});
            
            js_dir_error = rad2deg(circ_dist(deg2rad(js_dir{iState}), deg2rad(frme_vec{iState}.nominal_dir_deg)));
            joystick_acc = abs(1 - abs(js_dir_error) / 180);
            acc_detrend = joystick_acc - nanmean(joystick_acc);
            
            
            state_coh(cnt) = frme_vec{iState}.nominal_coh;
            
            [xc_acc_vec(cnt,:) lags] = xcorr(ecc_detrend(n_ofs:end),v_smooth_detrend(n_ofs:end),nLag,'normalized');
            [xc_ecc_vec(cnt,:) lags] = xcorr(acc_detrend(n_ofs:end),v_smooth_detrend(n_ofs:end),nLag,'normalized');
            [xc_ecc_acc(cnt,:) lags] = xcorr(ecc_detrend(n_ofs:end),acc_detrend(n_ofs:end),nLag,'normalized');
            [xc_control(cnt,:) lags] = xcorr(ecc_detrend(n_ofs:end-nshift+n_ofs),ecc_detrend(nshift:end),nLag,'normalized');
        end
    end
    
    snr_lst = unique(state_coh);
    if length(snr_lst)>7
        snr_lst = snr_lst(1:2:end);
    end
        
    for iCoh = 1:length(snr_lst)
        cidx = state_coh == snr_lst(iCoh);
        nanidx = ~isnan(sum(xc_ecc_acc,2));
        str{iCoh} = num2str(round(snr_lst(iCoh),2));
        avg_xc_acc_vec(iSubj,iCoh,:) = nanmean(xc_acc_vec(cidx' & nanidx,:));
        avg_xc_ecc_vec(iSubj,iCoh,:) = nanmean(xc_ecc_vec(cidx' & nanidx,:));
        avg_xc_ecc_acc(iSubj,iCoh,:) = nanmean(xc_ecc_acc(cidx' & nanidx,:));
        avg_xc_control(iSubj,iCoh,:) = nanmean(xc_control(cidx' & nanidx,:));
        
        clear tmp_corr
        tmp_corr = xc_ecc_acc(cidx' & nanidx,:);
        for iBlock = 1:size(tmp_corr)
            dat = tmp_corr(iBlock,:);
            try
                pk_pos{iSubj,iCoh}(iBlock)  = find(max(dat) == dat);
            catch
                vals = find(max(dat) == dat);
                pk_pos{iSubj,iCoh}(iBlock)  = vals(1);
            end
        end
    end
end

avg_pk_pos = cellfun(@median,pk_pos);
std_pk_pos = cellfun(@std,pk_pos);

%% PLOT
f                	= figure('units','centimeters','position',[0 0 5 5]);
ax                  = gca;
dat                 = avg_xc_ecc_acc;
lb_fs               = 8;
lw                  = 1.5;
col                 = cool(length(snr_lst));

ln = line([nLag nLag],[-1 1], 'Color', 'k','LineStyle', ':', 'LineWidth', 1);

for iCoh = 1:length(snr_lst)
    [ax, pl] = xc_plot_averages(ax,squeeze(dat(:,iCoh,:)),col(iCoh,:));
end

ax.XLabel.String = 'Lag [ms]';
ax.YLabel.String = 'XCorr coef';
ax.FontSize = lb_fs;
ax.XTick = [150-60 150 150+60 150+120];

for iLab = 1:length(ax.XTickLabel)
    ax.XTickLabel{iLab} = num2str( round((str2num(ax.XTickLabel{iLab})-nLag) * (1000/120)));
end

axis tight
ax.YLim = [-.1 .3];

dest_dir            = '/Users/fschneider/Documents/GitHub/CPR/Publications/2024_perceptual_confidence/FIG_solo_behaviour/raw/';
print(f, [dest_dir '/xcorr_joy_pop'], '-r500', '-dsvg');
print(f, [dest_dir '/xcorr_joy_pop'], '-r500', '-dpng');

%%
frme_ms             = 1000/120;

f                   = figure('units','centimeters','position',[0 0 5 5]); hold on
ax = subplot(2,1,1);
vl                  = violinplot( ((avg_pk_pos-nLag) * frme_ms) ./1e3);
vl                  = improveViolin(vl,col);
ax                  = gca;
ax.YLabel.String    = 'Median lag [s]';
ax.XTick            = 1:7;
ax.XTickLabel       = round(snr_lst,2)*100;
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
ax.XTickLabel       = round(snr_lst,2)*100;
ax.XLabel.String    = 'Coherence [%]';
ax.FontSize      	= lb_fs;
ax.Box              = 'off';
axis tight

dest_dir            = '/Users/fschneider/Documents/GitHub/CPR/Publications/2024_perceptual_confidence/FIG_solo_behaviour/raw/';
print(f, [dest_dir '/xcorr_joy_peaks'], '-r500', '-dsvg');
print(f, [dest_dir '/xcorr_joy_peaks'], '-r500', '-dpng');

%%
f                	= figure('units','centimeters','position',[0 0 15 15]);
lb_fs               = 8;
lw                  = 1.5;
col                 = cool(length(snr_lst));

for iPlot = 1:2
    
    ax          	= subplot(2,2,iPlot);
    
    if iPlot == 1
            dat  	= avg_xc_acc_vec;
            tstr    = 'Accuracy';
    else
            dat  	= avg_xc_ecc_vec;
            tstr    = 'Tilt';
    end
    
    ln              = line([nLag nLag],[-1 1], 'Color', 'k','LineStyle', ':', 'LineWidth', 1);
    
    for iCoh = 1:length(snr_lst)
        [ax, pl]    = xc_plot_averages(ax,squeeze(dat(:,iCoh,:)),col(iCoh,:));
    end
    
    ax.Title.String = tstr;
    ax.XLabel.String = 'Lag [ms]';
    
    if iPlot == 1
        ax.YLabel.String = {'Cross-correlation'; 'coefficient [norm.]'};
    end
    
    ax.FontSize     = lb_fs;
    ax.XTick        = [150-60 150 150+60 150+120];
    
    for iLab = 1:length(ax.XTickLabel)
        ax.XTickLabel{iLab} = num2str( round((str2num(ax.XTickLabel{iLab})-nLag) * (1000/120)));
    end
    
    axis tight
    ax.YLim         = [-.05 .05];
end

dest_dir            = '/Users/fschneider/Desktop/revision/';
print(f, [dest_dir '/xcorr_vec'], '-r500', '-dsvg');
print(f, [dest_dir '/xcorr_vec'], '-r500', '-dpng');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ax, pl] = xc_plot_averages(ax,dat,col)
axes(ax); hold on

lw                  = 1.5;
alp                 = .2;

% Boostrap confidence intervals
nRep        	= 5000;
[CI,~]        	= bootci(nRep,{@mean,dat},'Alpha',0.05);

% Prepare filled area
vec           	= 1:length(CI);
x_spacing     	= [vec fliplr(vec)];
ci            	= [CI(1,:) fliplr(CI(2,:))];

% Plot peak line
pk           	= find(mean(dat) == max(mean(dat)));
ln            	= line([pk pk],[-.1 max(mean(dat))], 'Color', col,'LineStyle', ':', 'LineWidth', lw/1.5);
   
% Overlay confidence intervals
fl           	= fill(x_spacing,ci,col,'EdgeColor','none', 'FaceAlpha', alp);

% Plot mean curve
pl              = plot(nanmean(dat), 'LineWidth', lw, 'Color', col);

end

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

% addpath '/Users/fschneider/Documents/GitHub/CPR/Matlab/'
% addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab/
% cd /Users/fschneider/Desktop/tmp/
% fname = '20220301_anm_CPRagent_block1_fxs.mwk2';
% d = MW_readFile(fname, 'include', {'AGNT_direction','#stimDisplay'}, 'dotPositions');

clearvars -except d

idx             = [];
idx.dot         = d.event == 'STIM_RDP_dot_positions';
idx.rdp_dir     = d.event == 'STIM_RDP_direction';
idx.rdp_coh     = d.event == 'STIM_RDP_coherence';
idx.agnt_dir    = d.event == 'AGNT_direction'; % STIM_AGNT_arc_direction

dp              = d.value(idx.dot);
dp_ts           = d.time(idx.dot);
rdp_dir         = d.value(idx.rdp_dir);
rdp_ts          = d.time(idx.rdp_dir);
rdp_coh         = d.value(idx.rdp_coh);
rdp_coh_ts      = d.time(idx.rdp_coh);
agnt_ts_tmp     = d.time(idx.agnt_dir);
agnt_dir_tmp    = d.value(idx.agnt_dir);
agnt_dir        = mod(cell2mat(agnt_dir_tmp(4:end)),360);
agnt_ts         = agnt_ts_tmp(4:end);

xIdx         	= logical(mod([1:size(dp{1},2)],2));
xpos{1}         = dp{1}(xIdx);
ypos{1}         = dp{1}(~xIdx);
c               = 1;
good_frme       = 1;
n               = length(dp);

for iFrme = 2:n
    disp(['Frame: ' num2str(iFrme)])
    
    if sum(isnan(dp{iFrme}))
        continue
    end
    
    % Dot position
    c                   = c+1;
    xIdx                = logical(mod([1:size(dp{iFrme},2)],2));
    xpos                = dp{iFrme}(xIdx);
    ypos                = dp{iFrme}(~xIdx);
    xpos_last           = dp{good_frme}(xIdx);      	% Last frame: x-position
    ypos_last         	= dp{good_frme}(~xIdx);      	% Last frame: y-position
    dts(c)              = dp_ts(iFrme);
    scoh(c)             = cell2mat(rdp_coh(find(dp_ts(iFrme) > rdp_coh_ts, 1,'last')));
    sdir(c)             = mod(cell2mat(rdp_dir(find(dp_ts(iFrme) > rdp_ts, 1,'last'))),360);
    
    for iDot = 2:length(ypos)
        vs              = [xpos_last(iDot), ypos_last(iDot)]; % dot position last frame
        ve              = [xpos(iDot), ypos(iDot)]; % dot position this frame
        dist(iDot)      = pdist([vs;ve],'euclidean'); % distance/vector length between points
        
        delta               = ve - vs; % vector that points from x1 to x2
        dot_dir{c}(iDot)	= mod(atan2d(delta(1),delta(2)),360); % Alternative: rad2deg(cart2pol(delta(2), delta(1)))
        dt{c}(iDot,:)       = delta;   
    end
    
    dot_idx{c}              = dist >= median(dist)-.0001 & dist <= median(dist)+.0001;
    good_frme               = iFrme;
end

for i = 2:n
    mdf                     = mean(dt{i}(dot_idx{i},:));                    % Mean x/y of all dots that didn't jump
    resultant(i)            = mod(atan2d(mdf(1),mdf(2)),360);      % Resultant vector [deg]
    res_length(i)           = pdist([[0 0];mdf],'euclidean');
end

figure;hold on
s1                          = stairs((dts-dts(1))./1e6,mod(sdir,360),'Linewidth',2);
p1                          = plot((dts-dts(1))./1e6,resultant,'Linewidth',2);
p2                          = plot((agnt_ts-dts(1))./1e6,agnt_dir,'Linewidth',2);
lg                          = legend([s1,p1,p2],'Veridical','Resultant','Agent') ;
lg.Location                 = 'south' ;
ax                          = gca;
ax.YLabel.String            = 'Direction [deg]';
ax.XLabel.String            = 'Time [s]';
% ax.XLim                     = [2460 2475];
ax.XLim                     = [2668 2685]; % CHECK COHERENCE 
ax.Title.String             = 'Example timecourse';
ax.Title.Interpreter        = 'none';
ax.FontSize                 = 16;

snr                         = unique(scoh);
lb                          = [];
dat                         = [];
for k = 1:length(snr)
    cidx                    = scoh == snr(k); 

    sig                     = abs(rad2deg(circ_dist(deg2rad(resultant(cidx)),deg2rad(sdir(cidx))))) < 1;
    r_coh(k)                = sum(sig)/length(sig);

    dat                     = [dat rad2deg(circ_dist(deg2rad(resultant(cidx)),deg2rad(sdir(cidx))))];
    lb                      = [lb repmat(k,1,sum(cidx))];
    rlen(k)                 = mean(res_length(cidx));
    avg(k)                  = mean(rad2deg(circ_dist(deg2rad(resultant(cidx)),deg2rad(sdir(cidx)))));
    sd(k)                   = std(rad2deg(circ_dist(deg2rad(resultant(cidx)),deg2rad(sdir(cidx)))));
end

f                           = figure;hold on
bx                          = boxplot(dat,lb, 'Color', [0 0 0]);
ax                          = gca;
ax.YLabel.String            = 'Deviation from veridical direction [deg]';
ax.XLabel.String            = 'Coherence [%]';
ax.Title.String             = ['Subject: AnM_agnt_block1: ' num2str(length(dat)) ' Frames'];
ax.Title.Interpreter        = 'none';
ax.XTick                    = 1:length(snr);
ax.XTickLabel               = snr;
ax.FontSize                 = 16;
print(f, '/Users/fschneider/Desktop/dir_deviation', '-r300', '-dpng');

% Check actual coherence level based on resultant vector
f                           = figure; hold on
pl                          = plot(r_coh, 'LineWidth',2,'Color','k');
ax                          = gca;
ax.YLabel.String            = '[%]';
ax.XLabel.String            = 'Coherence [%]';
ax.Title.String             = 'Diff(Resultant,Veridical) < 1deg';
ax.XTick                    = 1:length(snr);
ax.XTickLabel               = snr;
ax.FontSize                 = 16;
print(f, '/Users/fschneider/Desktop/resultant_coh', '-r300', '-dpng');

% Check coherence level based on dot direction
ofs = .01;
for m = 2:length(dot_dir)
    dots = dot_dir{m}(dot_idx{m}); 
    signal = abs(rad2deg(circ_dist(deg2rad(dots),deg2rad(sdir(m))))) < ofs;
    dcoh(m) = sum(signal)/length(signal);
end

f = figure; hold on
h1 = histogram(scoh,200,'FaceColor','k','FaceAlpha',.75);
h2 = histogram(dcoh,200,'FaceColor','r','FaceAlpha',.75);
ax                          = gca;
ax.YLabel.String            = 'No. Frame Pairs';
ax.XLabel.String            = 'Coherence [%]';
ax.Title.String             = 'Estimation of signal coherence based dot movements ';
ax.FontSize                 = 16;
print(f, '/Users/fschneider/Desktop/dots_coh', '-r300', '-dpng');


% Check evidence strength based on resultant vector length
f                           = figure; hold on
pl                          = plot(rlen, 'LineWidth',2,'Color','k');
ax                          = gca;
ax.YLabel.String            = 'Avg length of resultant vector';
ax.XLabel.String            = 'Coherence [%]';
ax.Title.String             = 'Evidence strength';
ax.XTick                    = 1:length(snr);
ax.XTickLabel               = snr;
ax.FontSize                 = 16;
print(f, '/Users/fschneider/Desktop/resultant_strength', '-r300', '-dpng');


for k = [5000 10000 15000 20000]
    figure; hold on
    histogram(mod(dot_dir{k},360),100)
    line([sdir(k) sdir(k)],[0 1000])
    title([num2str(sdir(k)) '___' num2str(scoh(k))])
end
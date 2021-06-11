close all

addpath /Users/fschneider/Documents/MATLAB/CircStat2012a
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab

load('/Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_psychophysics/Pilot_eye_fixation/mac_20210408_tbl.mat')

iTrl                = 1;                                                    % Trial
ex                  = t.trl_eye_x{iTrl};                                    % X data [dva]
ey                  = t.trl_eye_y{iTrl};                                    % Y data [dva]
ets                 = t.trl_eye_ts{iTrl};                                   % Timestamps
fx                  = t.trl_fix{iTrl};                                      % Fixation break flag
fx_ts               = t.trl_fix_ts{iTrl};                                   % Fixation break timestamps
xIdx                = ex < 4 & ex > -4;                                     % Index: Gaze outside of window
yIdx                = ey < 4 & ey > -4;
rdp_dir             = t.trl_rdp_dir{iTrl};                                  % RDP direction
rdp_dir_ts          = t.trl_rdp_dir_ts{iTrl};                               
rdp_coh             = t.trl_rdp_coh{iTrl};                                  % RDP coherence
rdp_coh_ts          = t.trl_rdp_coh_ts{iTrl};                                  

% Plot heat map
x                	= ex(xIdx & yIdx);
y                   = ey(xIdx & yIdx);
bins                = -4:.1:4;
f                   = figure; hold on

hist3([x',y'],'CDataMode','auto','FaceColor','interp','Ctrs',{bins bins});
view(2); axis equal; caxis([0 300]); colormap([0 0 0; jet(256)])

ax                  = gca;
ax.YLim             = [-4 4];
ax.XLim             = [-4 4];
ax.XLabel.String    = 'Horizontal location [dva]';
ax.YLabel.String    = 'Vertical location [dva]';
ax.FontSize         = 16;
cl                	= colorbar;

r                   = 4;
th                  = 0:pi/50:2*pi;
xunit               = r * cos(th);
yunit               = r * sin(th);
pl                  = plot(xunit, yunit, 'LineWidth',4,'Color', [1 1 1]);

% Plot polar distribution
f                   = figure;
[th,r]              = cart2pol(x,y);
pol                 = polarplot(th,r);
pol.Marker          = '.';
pol.LineStyle       = 'none';

%% Gaze accuracy
trlIdx              = cell2mat(cellfun(@length,t.trl_eye_ts,'UniformOutput',false)) > 1;
tt                  = t(trlIdx,:);
cc                  = 0;

for iTrl = 1:max(unique(t.trl_no))
    
ex                  = tt.trl_eye_x{iTrl};                                    % X data [dva]
ey                  = tt.trl_eye_y{iTrl};                                    % Y data [dva]
ets                 = tt.trl_eye_ts{iTrl};                                   % Timestamps
rdp_dir             = tt.trl_rdp_dir{iTrl};                                  % RDP direction
rdp_dir_ts          = tt.trl_rdp_dir_ts{iTrl};                               
[th,r]              = cart2pol(ex,ey);

for iState = 2:length(rdp_dir)
    if iState < length(rdp_dir)
        ssIdx       = ets > rdp_dir_ts(iState) & ets < rdp_dir_ts(iState+1);
    else
        ssIdx       = ets > rdp_dir_ts(iState);
    end
    
    clear x y excl gz_dir acc
    
    gz_dir          = mod(rad2deg(th(ssIdx)),360);
    
    excl            = r(ssIdx)>4;                                         	% Exclude fixation breaks [no time criterion here]
    gz_dir(excl)    = nan;
        
    df              = rad2deg(circ_dist(deg2rad(gz_dir),deg2rad(rdp_dir(iState))));	% Angular difference
    acc             = 1 - abs(df/180);                                              % Gaze accuracy
    
    cc = cc+1;
    avg_acc(cc)    	= nanmean(acc(end-149:end));                          	% Avg gaze accuracy in final 300ms of state (150 sample @ 500Hz) 
    
end
end

f                   = figure;
vs                  = violinplot(avg_acc,t.ss_coh);
ax                  = gca;
ax.YLabel.String    = 'Gaze accuracy [norm]';
ax.XLabel.String    = 'Coherence level';
ax.FontSize         = 16;

%% Model testing [GLM, LME]
clear ID coh rdp_dir js_dir acc fix

frv                 = load('/Users/fschneider/ownCloud/CPR_data/Pilot_free_viewing/CPR_pop_tbl.mat');
fix                 = load('/Users/fschneider/ownCloud/CPR_data/Pilot_eye_fixation/CPR_pop_tbl.mat');
t                   = [frv.t;fix.t];
fIdx                = logical([zeros(size(frv.t,1),1); ones(size(fix.t,1),1)]);
ID                  = t.ID;
coh                 = t.ss_coh;
rdp_dir             = t.rdp_dir;
js_dir              = cellfun(@(x) x(end-29:end),t.js_dir,'UniformOutput',false);   % Time window of interest [final 300ms, 30 samples @ 100Hz]
js_str              = cellfun(@(x) x(end-29:end),t.js_str,'UniformOutput',false);   % Time window of interest [final 300ms, 30 samples @ 100Hz]

for i = 1:length(js_dir)
    df              = rad2deg(circ_dist(deg2rad(js_dir{i}),deg2rad(rdp_dir(i))));   % Angular difference
    acc(i)          = mean(1 - abs(df/180));                                        % Motion tracking accuracy
    ecc(i)          = mean(js_str{i});                                              % Avg joystick eccentricity
end

tbl                 = table(ID, coh, acc', ecc', fIdx, 'VariableNames',{'ID','coh','acc','ecc', 'fix'});
tbl(strcmp(tbl.ID, 'piy_20210407') | strcmp(tbl.ID, 'stm_20210407'),:) = [];

% Intercept plus random effect for each level of the grouping variable g1: 'y ~ 1 + (1 | g1)'   
lme0                = fitlme(tbl,'acc~1 + ( 1 | ID)');
% Random intercept model with a fixed slope multiplying x1: 'y ~ x1 + (1 | g1)'
lme1                = fitlme(tbl,'acc~coh + ( 1 | ID)'); 
% Random intercepts and slopes, with possible correlation between them: 'y ~ x1 + (x1 | g1)' 
lme2                = fitlme(tbl,'acc~coh + ( coh | ID)'); 
% Independent random intercepts and slopes: 'y ~ x1 + (1 | g1) + (-1 + x1 | g1)'
lme3                = fitlme(tbl,'acc~coh + (1 | ID) + (-1 + coh | ID)'); 
% Random intercept model with independent main effects for g1 and g2, plus
% an independent interaction effect: 'y ~ 1 + (1 | g1) + (1 | g2) + (1 | g1:g2)'
lme                 = fitlme(tbl,'acc~1 + (1 | ID) + ( 1 | acc:ecc)'); 

stats               = compare(lme1,lme2);
[~,~,lmeStat]       = fixedEffects(lme);


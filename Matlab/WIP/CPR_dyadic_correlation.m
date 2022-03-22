addpath /Users/fschneider/Documents/MATLAB/CircStat2012a
addpath /Users/fschneider/Documents/MATLAB/cbrewer

a = load('/Users/fschneider/Desktop/tmp8/block2/20220303_anm_cpr_tbl.mat');
b = load('/Users/fschneider/Desktop/tmp8/block2/20220303_cah_cpr_tbl.mat');

snr                     = unique(a.t.ss_coh);
nSamples                = 15;
cl                      = linspace(0,1,length(snr));
figure

% FILTER
% only hit states
% only right after change
% -> reduce variability
% FINDING: Differences between subjects
% What produces this type of effect? CLUSTER SUBJECTS! Confusion/Detrimental performance vs Integration of additional
% information - better
% Correelation between str and acc fro different classes
% look also at agent vs dyadic cluster according to psychometric curve -> also modulated by signal2noise
% follower vs leader vs independent - might change with SNR
% decision times for each subject

for iCoh = 1:length(snr)
    clear cIdx rdp_dir t.js_dir
    
    cIdx = b.t.ss_coh == snr(iCoh);
    cIdx(cellfun(@length,t.js_str) < 100) = false;
    
    % Joystick displacement
    str = b.t.js_str(cIdx);
    
    % Joystick accuracy
    rdp_dir                 = a.t.rdp_dir(cIdx);
    js_dir                  = a.t.js_dir(cIdx);
    
    for iState = 1:length(rdp_dir)
        clear js_dev
        js_dev              = rad2deg(circ_dist(deg2rad(js_dir{iState}),deg2rad(rdp_dir(iState))));  % Minimum RDP-Joystick difference
        js_acc{iCoh}(iState)      = nanmean(abs(1 - abs(js_dev(end-nSamples:end)) / 180));                           % Joystick accuracy
        js_str{iCoh}(iState)      = nanmean(str{iState}(end-nSamples:end));
    end  
    
    scatter(js_acc{iCoh},js_str{iCoh},'MarkerFaceColor', [cl(iCoh) 0 0],'MarkerEdgeColor', [cl(iCoh) 0 0])
    hold on
end


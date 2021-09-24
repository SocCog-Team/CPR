function [] = CPR_fit_sigmoid(cpr,afc)

addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a
addpath /Users/fschneider/Documents/MATLAB/palamedes1_10_9/Palamedes

%% CONTINUOUS DATA: Binary

sub                     = unique(cpr.ID);                                   	% Subject ID
snr                     = unique(cpr.ss_coh);                                 % Stimulus coherence level
snr(snr == 0)           = [];
nSamples                = 9;                                               	% Integration period [1 sample == 10ms]
tresh                   = 45;                                              	% +/-45deg --> 90 degree tolerance == 25% chance level

for iSS = 2:size(cpr.js_dir,1)
    if isnan(cpr.js_dir{iSS})
        js_dir(iSS,:)   = nan;
        continue
    else
        js_dir(iSS,:)	= median(cpr.js_dir{iSS}(end-nSamples:end));
    end
end

dff                     = rad2deg(circ_dist(deg2rad(cpr.rdp_dir),deg2rad(js_dir)));  % Raw RDP-Joystick difference
tdff                    = abs(dff) < tresh;

for iCoh = 1:length(snr)
    clear cIdx
    cIdx                = cpr.ss_coh == snr(iCoh);
    HitNo(1,iCoh)       = sum(tdff(cIdx));
    OutOfNum(1,iCoh)    = length(tdff(cIdx));
    hir(1,iCoh)         = HitNo(1,iCoh) / OutOfNum(1,iCoh);
end

%% CONTINUOUS DATA: Precision-based

sub                     = unique(cpr.ID);                                   	% Subject ID
snr                     = unique(cpr.ss_coh);                                 % Stimulus coherence level
snr(snr == 0)           = [];
nSamples                = 9;                                               	% Integration period [1 sample == 10ms]

for iSS = 2:size(cpr.js_dir,1)
    if isnan(cpr.js_dir{iSS})
        js_dir(iSS,:)   = nan;
        continue
    else
        js_dir(iSS,:)	= median(cpr.js_dir{iSS}(end-nSamples:end));
    end
end

dff                     = rad2deg(circ_dist(deg2rad(cpr.rdp_dir),deg2rad(js_dir)));  % Raw RDP-Joystick difference

for iCoh = 1:length(snr)
    clear cIdx
    cIdx                = cpr.ss_coh == snr(iCoh);
    mdiff(iCoh)         = nanmean(abs(dff(cIdx)));
end

mprc                    = 1-mdiff./180;

%% 4AFC DATA
afc = [afc;afc;afc;afc];
validTrls               = strcmp(afc.trl_outcome, 'wrong') | strcmp(afc.trl_outcome, 'hit');
HI                      = strcmp(afc.trl_outcome, 'hit');
snr_4afc                = unique(afc.trl_coh);

for iCoh = 1:length(snr_4afc)
    cidx                = afc.trl_coh == snr_4afc(iCoh);
    HitNo(2,iCoh)       = sum(HI(cidx));
    OutOfNum(2,iCoh)    = sum(validTrls(cidx));
    hir(2,iCoh)       	= HitNo(2,iCoh) / OutOfNum(2,iCoh);
end

%% FIT & PLOT

for iExp = 1:2
    
    clear vec
    if iExp == 1
        vec = snr';
        cl = [0 0 0];
        disp('CPR Fit: Binary')
    else
        vec = snr_4afc';
        cl = [1 0 0];
        disp('4AFC Fit')
    end
    
    %%% SEE PALAMEDES DEMO %%%
    % 'PF': psychometric function to be fitted. Passed as an inline function.
    % Options include:
    %	@PAL_Logistic
    % 	@PAL_Weibull
    %  	@PAL_Gumbel (i.e., log-Weibull)
    %  	@PAL_Quick
    %  	@PAL_logQuick
    % 	@PAL_CumulativeNormal
    %  	@PAL_HyperbolicSecant
    PF = @PAL_Logistic;
    
    %Threshold and Slope are free parameters, guess and lapse rate are fixed
    paramsFree = [0 1 0 0];  %1: free parameter, 0: fixed parameter
    
    %Parameter grid defining parameter space through which to perform a
    %brute-force search for values to be used as initial guesses in iterative
    %parameter search.
    searchGrid.alpha = vec(1):.001:vec(end); % threshold
    searchGrid.beta = logspace(0,2,101); % slope
    searchGrid.gamma = 0.25;  % guess-rate
    searchGrid.lambda = 0.02;  % lapse-rate
    
    %Perform fit
    disp('Fitting function.....');
    [paramsValues LL exitflag] = PAL_PFML_Fit(vec,HitNo(iExp,1:length(vec)), ...
        OutOfNum(iExp,1:length(vec)),searchGrid,paramsFree,PF);
    
    disp('done:')
    message = sprintf('Threshold estimate: %6.4f',paramsValues(1));
    disp(message);
    message = sprintf('Slope estimate: %6.4f\r',paramsValues(2));
    disp(message);
    
    %Number of simulations to perform to determine standard error
    B=400;
    disp('Determining standard errors.....');
    
    [SD paramsSim LLSim converged] = PAL_PFML_BootstrapNonParametric(...
        vec,HitNo(iExp,1:length(vec)), OutOfNum(iExp,1:length(vec)), [], paramsFree, B, PF,...
        'searchGrid',searchGrid);
    
    disp('done:');
    message = sprintf('Standard error of Threshold: %6.4f',SD(1));
    disp(message);
    message = sprintf('Standard error of Slope: %6.4f\r',SD(2));
    disp(message);
    
    %Number of simulations to perform to determine Goodness-of-Fit
    B=1000;
    disp('Determining Goodness-of-fit.....');
    
    [Dev pDev] = PAL_PFML_GoodnessOfFit(vec,HitNo(iExp,1:length(vec)), OutOfNum(iExp,1:length(vec)), ...
        paramsValues, paramsFree, B, PF, 'searchGrid', searchGrid);
    
    disp('done:');
    message = sprintf('Deviance: %6.4f',Dev);
    disp(message);
    message = sprintf('p-value: %6.4f',pDev);
    disp(message);
    
    %Create simple plot
    vec = round(vec,2);
    ProportionCorrectObserved=HitNo(iExp,1:length(vec))./OutOfNum(iExp,1:length(vec));
%     StimLevelsFineGrain=[min(vec):max(vec)./1000:max(vec)];
    StimLevelsFineGrain=[0:.01:1];
    ProportionCorrectModel = PF(paramsValues,StimLevelsFineGrain);
    
    plot(StimLevelsFineGrain,ProportionCorrectModel,'-','color',cl,'linewidth',4);
    plot(vec,ProportionCorrectObserved,'color',cl,'marker', '.','LineStyle','none','markersize',40);
    set(gca, 'fontsize',16);
    set(gca, 'Xtick',[0:.2:.8]);
    xlabel('Coherence');
    ylabel('% correct'); 
    ylim([0 1])
    xlim([0 .8])
end

% %%% Precision-based addition to plot %%%
% yyaxis right
% ft                  = fit(snr,mprc','1./(1+exp(b1+b2*x))');
% rnge                = 0:.01:.8;
% p                   = plot(rnge,ft(rnge), 'LineStyle', '-','linewidth',4, 'color',[.5 .5 .5]);
% sc                  = scatter(snr,mprc);
% sc.MarkerFaceColor  = [.5 .5 .5];
% sc.MarkerEdgeColor  = [.5 .5 .5];
% ax                  = gca;
% ax.YAxis(1).Color   = [0 0 0];
% ax.YAxis(2).Color   = [.5 .5 .5];
% ax.YLabel.String    = 'Precision'; 

end
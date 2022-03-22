function f = CPR_fit_sigmoid(snr_cpr,snr_4afc, HitNo, OutOfNum, macc)

addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a
addpath /Users/fschneider/Documents/MATLAB/palamedes1_10_9/Palamedes

%% FIT & PLOT

f                   = figure; hold on
lw                  = 4;
alph                = .7;
cl_acc              = [.85 0 0];
cl_perc             = [0 .6 .8];

%%% Precision-based addition to plot %%%
yyaxis right
ft                  = fit(snr_cpr,macc,'1./(1+exp(b1+b2*x))');
rnge                = 0:.01:.8;
p                   = plot(rnge,ft(rnge), 'LineStyle', '-','linewidth',lw, 'color',[cl_acc alph]);
sc                  = scatter(snr_cpr,macc);
sc.MarkerFaceColor  = cl_acc;
sc.MarkerEdgeColor  = cl_acc;
sc.SizeData         = 100;
ax                  = gca;
ax.YAxis(1).Color   = cl_perc/2;
ax.YAxis(2).Color   = cl_acc;
ax.YLabel.String    = 'Accuracy'; 
ax.YLim             = [0 1]; 

yyaxis left

for iExp = 1:2
    
    clear vec
    if iExp == 1
        vec = snr_cpr';
        cl = [cl_perc alph];
        disp('CPR Fit: Binary')
    else
        vec = snr_4afc';
        cl = [cl_perc/2 alph];
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
    
    pl(iExp) = plot(StimLevelsFineGrain,ProportionCorrectModel,'-','color',cl,'linewidth',lw);
    pd(iExp) = plot(vec,ProportionCorrectObserved,'color',cl,'marker', '.','LineStyle','none','markersize',40);
    set(gca, 'fontsize',16);
    set(gca, 'Xtick',[0:.2:.8]);
    xlabel('Coherence');
    ylabel('% correct'); 
    ylim([0 1])
    xlim([0 .8])
end

ln                  = line([0 .8],[.25 .25]);
ln.LineStyle        = '--';
ln.LineWidth        = 2;
ln.Color            = [0 0 0];

lg                  = legend([p,pl(1),pl(2),ln], 'CPR_acc','CPR_thresh','4AFC','Chance');
lg.Location         = 'southeast';
lg.Interpreter      = 'none';
end
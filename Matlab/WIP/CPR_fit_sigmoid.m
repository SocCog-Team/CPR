function f = CPR_fit_sigmoid(snr, HitNo, OutOfNum)

addpath /Users/fschneider/Documents/MATLAB/palamedes1_10_9/Palamedes

%% FIT & PLOT

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
searchGrid.alpha = snr(1):.001:snr(end); % threshold
searchGrid.beta = logspace(0,2,101); % slope
searchGrid.gamma = 0.25;  % guess-rate
searchGrid.lambda = 0.02;  % lapse-rate

%Perform fit
disp('Fitting function.....');
[paramsValues LL exitflag] = PAL_PFML_Fit(snr,HitNo,OutOfNum,searchGrid,paramsFree,PF);

disp('done:')
message = sprintf('Threshold estimate: %6.4f',paramsValues(1));
disp(message);
message = sprintf('Slope estimate: %6.4f\r',paramsValues(2));
disp(message);

%Number of simulations to perform to determine standard error
B=400;
disp('Determining standard errors.....');

[SD paramsSim LLSim converged] = PAL_PFML_BootstrapNonParametric(snr,HitNo, OutOfNum, [], paramsFree, B, PF,'searchGrid',searchGrid);

disp('done:');
message = sprintf('Standard error of Threshold: %6.4f',SD(1));
disp(message);
message = sprintf('Standard error of Slope: %6.4f\r',SD(2));
disp(message);

%Number of simulations to perform to determine Goodness-of-Fit
B=1000;
disp('Determining Goodness-of-fit.....');

[Dev pDev] = PAL_PFML_GoodnessOfFit(snr,HitNo,OutOfNum,paramsValues, paramsFree, B, PF, 'searchGrid', searchGrid);

disp('done:');
message = sprintf('Deviance: %6.4f',Dev);
disp(message);
message = sprintf('p-value: %6.4f',pDev);
disp(message);

%Create simple plot
snr = round(snr,2);
ProportionCorrectObserved=HitNo./OutOfNum;
%     StimLevelsFineGrain=[min(vec):max(vec)./1000:max(vec)];
StimLevelsFineGrain=[0:.01:1];
ProportionCorrectModel = PF(paramsValues,StimLevelsFineGrain);

%% PLOT

f                   = figure; hold on
lw                  = 4;
alph                = 1;

pl                  = plot(StimLevelsFineGrain,ProportionCorrectModel,'-','color',[0 0 0 alph],'linewidth',lw);
pd                  = plot(snr,ProportionCorrectObserved,'color',[0 0 0 alph],'marker', '.','LineStyle','none','markersize',40);
ax                  = gca;
ax.FontSize         = 16;
ax.XTick            = [0:.2:.8];
ax.XLabel.String    = 'Coherence';
ax.YLabel.String    = '% correct';
ax.YLim             = [0 1];
ax.XLim             = [0 .8];
ax.Title.String     = '4AFC psychometric curve';
ln                  = line(ax.XLim,[.25 .25],'Color', 'k', 'LineStyle', '--');

end

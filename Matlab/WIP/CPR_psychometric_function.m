function [snr, dat, model] = CPR_psychometric_function(pth, fname, dest_dir)

% Import 4AFC data
var_import = {
    'RDP_', ...
    'TRIAL_'};

d = MW_readFile([pth fname], 'include', var_import,'~typeOutcomeCheck'); 	% Import .mwk2 sesion file
idx.tOn                 = d.event == 'TRIAL_start';
idx.tEnd                = d.event == 'TRIAL_end';
idx.RDP_coh             = d.event == 'RDP_coherence';
idx.outcome             = d.event == 'TRIAL_outcome';

% Trial timestamps
trl.tOn                 = d.time(idx.tOn);
trl.tEnd                = d.time(idx.tEnd);

% Remove trials with duration <30us
excl                    = trl.tEnd - trl.tOn < 30;
trl.tOn(excl)           = [];
trl.tEnd(excl)          = [];

for iTrl = 1:length(trl.tEnd)
    clear trlIdx
    trlIdx            	= d.time >= trl.tOn(iTrl) & d.time <= trl.tEnd(iTrl);
    trl.coh(iTrl)     	= getTrialData(d.value, trlIdx, idx.RDP_coh);
    trl.outcome{iTrl} 	= getTrialData(d.value, trlIdx, idx.outcome);
end

%% Extract 4AFC data

ER                      = cell2mat(cellfun(@(x)strcmp(x,'wrong'),trl.outcome,'UniformOutput', false));
HI                      = cell2mat(cellfun(@(x)strcmp(x,'hit'),trl.outcome,'UniformOutput', false));
validTrls               = ER | HI;
snr                     = unique(trl.coh);
snr(isnan(snr))         = [];

for iCoh = 1:length(snr)
    cidx                = trl.coh == snr(iCoh);
    HitNo(iCoh)         = sum(HI(cidx));
    OutOfNum(iCoh)      = sum(validTrls(cidx));
    hir(iCoh)       	= HitNo(iCoh) / OutOfNum(iCoh);
end

%% Fit and plot psychometric function

[f, dat, model]         = CPR_fit_sigmoid(snr, HitNo, OutOfNum);

print(f, [dest_dir fname(1:16)], '-r300', '-dpng');

end
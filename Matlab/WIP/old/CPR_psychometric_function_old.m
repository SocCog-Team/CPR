function CPR_psychometric_function(pth,fname_4afc,fname_cpr_tbl)

cd([pth 'raw/'])

% Import 4AFC data
var_import = {
    'RDP_', ...
    'TRIAL_'};

d = MW_readFile([fname_4afc], 'include', var_import,'~typeOutcomeCheck'); 	% Import .mwk2 sesion file
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

%% Extract CPR data

load([pth 'summary/' fname_cpr_tbl])

sub                     = unique(t.ID);                                 	% Subject ID
snr_cpr                	= unique(t.ss_coh);                                 % Stimulus coherence level
nSamples                = 14;                                               % Integration period [1 sample == 10ms]
tresh                   = 45;                                              	% +/-45deg --> 90 degree tolerance == 25% chance level

for iSS = 1:size(t.js_dir,1)
    if isnan(t.js_dir{iSS})
        js_dir(iSS,:)   = nan;
        continue
    else
        js_dir(iSS,:)	= median(t.js_dir{iSS}(end-nSamples:end));
    end
end

dff                     = rad2deg(circ_dist(deg2rad(t.rdp_dir),deg2rad(js_dir)));  % Raw RDP-Joystick difference
tdff                    = abs(dff) < tresh;

% Threshold-based
for iCoh = 1:length(snr_cpr)
    clear cIdx
    cIdx                = t.ss_coh == snr_cpr(iCoh);
    HitNo(1,iCoh)       = sum(tdff(cIdx));
    OutOfNum(1,iCoh)	= length(tdff(cIdx));
    hir(1,iCoh)         = HitNo(1,iCoh) / OutOfNum(1,iCoh);
end

% Accuracy-based
for iCoh = 1:length(snr_cpr)
    clear cIdx
    cIdx                = t.ss_coh == snr_cpr(iCoh);
    mdiff(iCoh)         = nanmean(abs(dff(cIdx)));
end
macc                    = 1-mdiff./180;

%% Extract 4AFC data

ER                      = cell2mat(cellfun(@(x)strcmp(x,'wrong'),trl.outcome,'UniformOutput', false));
HI                      = cell2mat(cellfun(@(x)strcmp(x,'hit'),trl.outcome,'UniformOutput', false));
validTrls               = ER | HI;
snr_4afc                = unique(trl.coh);
snr_4afc(isnan(snr_4afc)) = [];

for iCoh = 1:length(snr_4afc)
    cidx                = trl.coh == snr_4afc(iCoh);
    HitNo(2,iCoh)       = sum(HI(cidx));
    OutOfNum(2,iCoh)    = sum(validTrls(cidx));
    hir(2,iCoh)       	= HitNo(2,iCoh) / OutOfNum(2,iCoh);
end

%% Fit psychometric function

f                       = CPR_fit_sigmoid(snr_cpr,snr_4afc', HitNo, OutOfNum, macc');
ax                      = gca;
ax.Title.String         = fname_4afc(10:12);

print(f, [pth 'summary/summary_' fname_4afc(10:12) '_psychometric_curves'], '-r300', '-dpng');

end
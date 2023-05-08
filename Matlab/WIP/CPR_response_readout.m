% Target score
function out = CPR_response_readout(in, nSample)

c                           = 0;
snr                         = unique(in.rdp_coh);

% Target response params
for iState = 1:size(in.trg_ts,1)
    for iTrg = 1:length(in.trg_ts{iState})
        c                   = c +1;
        score_cum(c)        = in.trg_score{iState}(iTrg);
        score_coh(c)        = in.rdp_coh(iState);
        score_hi(c)         = in.trg_hit{iState}(iTrg);
    end
end

tmp_score                   = score_cum(~isnan(score_cum));
tscore(~isnan(score_cum))   = [0 diff(tmp_score)];
tscore(tscore < 0)          = nan;

for iCoh = 1:length(snr)
    clear cIdx rdp_dir js_dir
    
    % Coherence index
    cIdx = in.rdp_coh == snr(iCoh);
    cIdx(cellfun(@length,in.js_str) < 100) = false;
    
    % Hit rate
    nhi                     = sum(cellfun(@sum,in.trg_hit(cIdx)));
    ntrg                    = sum(cellfun(@numel,in.trg_hit(cIdx)));
    out.hir(iCoh)           = nhi / ntrg;
    
    % Target score
    out.trg_mscore(iCoh)	= nanmean(tscore(score_coh  == snr(iCoh) & score_hi == true));
    out.trg_score{iCoh}  	= tscore(score_coh  == snr(iCoh) & score_hi == true);
    
    % Joystick displacement
    out.mstr(iCoh)         	= nanmedian(cellfun(@(x) nanmedian(x(end-nSample:end)), in.js_str(cIdx)));
    out.str{iCoh}           = cellfun(@(x) nanmedian(x(end-nSample:end)), in.js_str(cIdx));
    
    % Joystick accuracy before first target
    t1_ts                   = cellfun(@(x) x(1), in.trg_ts);
    f1_ts                   = cellfun(@(x) x(1), in.frme_ts);
    trgIdx                  = (t1_ts-f1_ts) > 1e6;
    rdp_dir                 = in.rdp_dir(cIdx & in.trg_shown & trgIdx);
    js_dir                  = in.js_dir(cIdx & in.trg_shown & trgIdx);
    frmes                   = in.frme_ts(cIdx & in.trg_shown & trgIdx);
    trg1_ts                 = t1_ts(cIdx & in.trg_shown & trgIdx);
    
    for iState = 1:length(rdp_dir)
        clear js_dev
        smpl_idx            = find(frmes{iState} < trg1_ts(iState),1,'last')-nSample : find(frmes{iState} < trg1_ts(iState),1,'last');
        js_dev              = rad2deg(circ_dist(deg2rad(js_dir{iState}(smpl_idx)),deg2rad(rdp_dir(iState))));  % Minimum RDP-Joystick difference
        js_acc(iState)      = nanmean(abs(1 - abs(js_dev) / 180));                           % Joystick accuracy
    end
    
    out.carr(iCoh)         	= snr(iCoh);
    out.macc(iCoh)        	= nanmedian(js_acc);
    out.acc{iCoh}        	= js_acc;
    
end
end

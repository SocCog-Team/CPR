% Add relevant directories
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/WIP/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/Helper_functions/
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab
addpath /Users/fschneider/Documents/MATLAB/cbrewer/

close all
clear all

pth                                     = '/Volumes/T7_Shield/CPR_psychophysics/';
fname                                   = 'Subjects_summary.xlsx';
x                                       = readtable([pth fname]);
sbj_lst                                 = x.Abbreviation;
dte_format                              = 'yyyymmdd';
sbj_lst(cellfun(@isempty,sbj_lst))      = [];
nLag                                    = 150;
nSample                                 = 30;
[solo_score, agnt_score, dyad_score]    = getFinalScore(pth, sbj_lst);

sbj_lst(cellfun(@(x) strcmp(x, 'RoH'),sbj_lst)) = [];

for iSubj = 1:length(sbj_lst)
    
    disp(['Processing subject: ' sbj_lst{iSubj}])
    
    % Import reaction time
    load([pth sbj_lst{iSubj} '/summary/' sbj_lst{iSubj} '_RT.mat'])
    mrt(iSubj)                          = mean(rt.dat);
    id_solo{iSubj}                    	= sbj_lst{iSubj};
    
    load([pth sbj_lst{iSubj} '/summary/' sbj_lst{iSubj} '_SNRfit.mat'])
    psy_df                              = abs(psy_func.model - .5);
    tresh_smpl(iSubj)                 	= find(psy_df == min(psy_df));
    thresh_snr(iSubj)                   = psy_func.model_snr(tresh_smpl(iSubj));
    
    % Import solo performance
    data_pth                            = [pth sbj_lst{iSubj} '/summary/'];
    solo_tbl                            = [];
    agnt_tbl                            = [];
    tmp_solo                            = [];
    tmp_agnt                           	= [];
    cc                                  = 0;
    ccc                                 = 0;
    
    if isdir(data_pth)
        cd(data_pth)
        mat_files                       = dir('*.mat');
        
        for iFile = 1:length(mat_files)
            if contains(mat_files(iFile).name,'CPRsolo')
                tmp_tbl                 = load(mat_files(iFile).name);
                
                % Select sessions with proper coherence level
                if sum(unique(tmp_tbl.t.rdp_coh) < .1) > 2
                    continue
                end
                
                % Concatenate session summary tables
                solo_tbl                    = [solo_tbl; tmp_tbl.t];
            end
            
            if contains(mat_files(iFile).name,'CPRagent') && contains(mat_files(iFile).name,lower(sbj_lst{iSubj}))
                tmp_tbl                     = load(mat_files(iFile).name);
                
                % Concatenate session summary tables
                agnt_tbl                    = [agnt_tbl; tmp_tbl.t];
            end
        end
        
        if ~isempty(solo_tbl)
            % Performance analysis
            solo_perf(iSubj)                = response_readout(solo_tbl, nSample);
            
            % Extract trials for correlation analysis
            trl_break                       = [0; find(diff(solo_tbl.trl_no) ~= 0)];
            [tmp_solo, tc]                  = extractTrials(trl_break, solo_tbl, tmp_solo, cc);
            
            % Correlation analysis
            [cr_solo{iSubj},ps]             = CPR_correlation_analysis_WIP(tmp_solo, nLag, false);
        end
        
        if ~isempty(agnt_tbl)
            % Performance analysis
            agnt_perf(iSubj)                = response_readout(agnt_tbl, nSample);
            
            % Extract trials for correlation analysis
            trl_break                       = [0; find(diff(agnt_tbl.trl_no) ~= 0)];
            [tmp_agnt, tc]                  = extractTrials(trl_break, agnt_tbl, tmp_agnt, ccc);
            
            % Correlation analysis
            [cr_agnt{iSubj},ps]             = CPR_correlation_analysis_WIP(tmp_agnt, nLag, false);
        end
    end
end

%% Dyadic sessions

c = 0;

for iDyad = 19:54
    
    disp(['Processing dyad: ' num2str(iDyad)])
    
    % Specify directories
    summ_pth                            = [pth 'Dyad' num2str(iDyad) '/summary/'];
    
    % Extract .mwk2 file names from source directory
    cd(summ_pth)
    mat_files                           = dir('*.mat');
    
    if isempty(mat_files)
        continue
    end
    
    c                                   = c+1;
    p1                                  = [];
    p2                                  = [];
    
    for iFile = 1:length(mat_files)
        tbl                             = [];
        tbl                             = load(mat_files(iFile).name);
        
        if iFile == 1 || iFile == 2
            p1                          = [p1; tbl.t];
        else
            p2                          = [p2; tbl.t];
        end
    end
    
    if ~isempty(p1)
        % Performance analysis
        dyad_perf{c,1}                  = response_readout(p1, nSample);
        dyad_perf{c,2}              	= response_readout(p2, nSample);
        
        % Extract trials for correlation analysis
        trl_break                       = [0; find(diff(p1.trl_no) ~= 0)];
        
        for iSubj = 1:2
            
            if iSubj == 1
                dat_tbl                 = p1;
            else
                dat_tbl                 = p2;
            end
            
            cc                          = 0;
            tmp                         = {};
            [tmp, cc]                   = extractTrials(trl_break, dat_tbl, tmp, cc);
            
            % Correlation analysis
            [cr_dyad{c,iSubj},ps]       = CPR_correlation_analysis_WIP(tmp, nLag, false);
            id_dyad{c}                  = [p1.ID(1) p2.ID(1)];
        end
    end
end


%% PLOT

f                           = figure('units','normalized','position',[0 0 .5 1]);
clm                         = linspace(.1,.77,3);
height                      = [.79 .5 .24 .08];
clmns                      	= linspace(.09,.8,4);
lb_fs                       = 14;
lg_fs                       = 10;
lw                          = 3;
dim                         = [.2 .2];
alp                         = .35;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Solo vs Dyadic score %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax0                       	= axes('Position', [clm(1) height(1) dim]); hold on
ln                          = line([100 300], [100 300]);
ln.LineStyle                = '--';
ln.LineWidth                = lw/2;
ln.Color                    = [0 0 0];

for i = 1:length(agnt_score)
    sc_agnt(i)            	= scatter(mean(agnt_score{i}),mean(solo_score{i}), 'MarkerFaceColor', [.5 .5 .5],'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .75);
    sc_dyad(i)            	= scatter(mean(dyad_score{i}),mean(solo_score{i}), 'MarkerFaceColor', [0 0 0],'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .75);
end

ax0.XLabel.String           = {'Avg reward score', '[Dyad]'};
ax0.YLabel.String           = {'Avg reward score', '[Solo]'};
ax0.XLim                    = [100 300];
ax0.YLim                    = [100 300];
ax0.FontSize                = lb_fs;

lg                          = legend([sc_agnt(1) sc_dyad(1)], 'AGNT','HUMAN' );
lg.Location                 = 'southeast';
lg.FontSize                 = lg_fs;
lg.Box                      = 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Hit rate comparison %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iSub = 1:length(solo_perf)
    hir.solo(iSub)          = solo_perf(iSub).hir_pool;
    if ~isempty(agnt_perf(iSub).hir_pool)
        hir.agnt(iSub)   	= agnt_perf(iSub).hir_pool;
    else
        hir.agnt(iSub) 	= nan;
  	end
end

hir.dyad = [];
for iDyad = 1:size(dyad_perf,1)
    hir.dyad                = [hir.dyad dyad_perf{iDyad,1}.hir_pool];
    hir.dyad                = [hir.dyad dyad_perf{iDyad,2}.hir_pool];
end

nRep = 500;
e.solo                      = bootci(nRep, {@mean, hir.solo},'alpha', .001);
e.agnt                      = bootci(nRep, {@nanmean, hir.agnt},'alpha', .001);
e.dyad                      = bootci(nRep, {@mean, hir.dyad},'alpha', .001);

ax1                      	= axes('Position', [clm(2) height(1) dim]); hold on
bp                          = bar([1:3],[mean(hir.solo) nanmean(hir.agnt) mean(hir.dyad)]);
bp.FaceColor                = [.5 .5 .5];
bp.EdgeColor                = 'none';

er                          = errorbar([1:3],[mean(hir.solo) nanmean(hir.agnt) mean(hir.dyad)],[mean(hir.solo)-e.solo(2) nanmean(hir.agnt)-e.agnt(2) mean(hir.dyad)-e.dyad(2)],[e.solo(1)-mean(hir.solo) e.agnt(1)-nanmean(hir.agnt) e.dyad(1)-mean(hir.dyad)]);    
er.Color                    = [0 0 0];                            
er.LineStyle                = 'none';
er.LineWidth                = lw/1.5;

ax1.YLabel.String           = 'Avg hit rate';
ax1.XLim                    = [.5 3.5];
ax1.YLim                    = [.5 .68];
ax1.YTick                 	= .5:.05:.7;
ax1.XTickLabel              = {'Solo','Agnt','Dyad'};
ax1.FontSize                = lb_fs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Timing - performance correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iD = 1:size(cr_dyad,1)
    %%% Timing leader
    lgs(iD,:)             	= [median(cr_dyad{iD,1}.lag) median(cr_dyad{iD,2}.lag)];
    [lgs(iD,:), lgs_lead(iD,:)]= sort(lgs(iD,:));
    
    %%% Performance leader
    if (nanmean(cell2mat(dyad_perf{iD,2}.trg_score)) - nanmean(cell2mat(dyad_perf{iD,1}.trg_score))) < 0
        score_lead(iD,:)   	= [1 2];
    elseif (nanmean(cell2mat(dyad_perf{iD,2}.trg_score)) - nanmean(cell2mat(dyad_perf{iD,1}.trg_score))) > 0
        score_lead(iD,:)  	= [2 1];
    end  
    
    for iS = 1:2
        try
            sIdx                = cellfun(@(x) strcmp(convertStringsToChars(id_dyad{iD}(lgs_lead(iD,iS))),x),cellfun(@lower,id_solo,'uni',false));
            afc4_tresh(iD,iS)   = thresh_snr(sIdx);
            avg_acc(iD,iS)      = nanmean(cell2mat(solo_perf(sIdx).acc));
        catch
            afc4_tresh(iD,iS)   = nan;
            avg_acc(iD,iS)      = nan;
        end
    end
end

% Sort - there must be a single-line option
for k = 1:size(avg_acc,1)
tmp_acc(k,:) = avg_acc(k,lgs_lead(k,:));
end

ax5                      	= axes('Position', [clm(3) height(1) dim]); hold on
ex_idx                      = logical(sum(isnan(tmp_acc),2));
xdat                        = tmp_acc(~ex_idx,:);               
ydat                        = lgs(~ex_idx,:); % already sorted
dcol                        = {[0 0 0], [.5 .5 .5]};
lcol                        = jet(size(xdat,1));
alph                        = .8;

for iD = 1:size(xdat,1)   
    ln                      = line([xdat(iD,1) xdat(iD,2)],[ydat(iD,1) ydat(iD,2)]);
    ln.LineWidth            = 2;
    ln.Color             	= [lcol(iD,:) alp];
    
    sc                   	= scatter(xdat(iD,1), ydat(iD,1));
    sc.MarkerEdgeColor    	= 'none';
    sc.MarkerFaceColor     	= [dcol{1}];
    sc.MarkerFaceAlpha     	= alph;
    sc.MarkerFaceAlpha     	= alph;
    sc.MarkerFaceAlpha     	= alph;

    sc                     	= scatter(xdat(iD,2), ydat(iD,2));
    sc.MarkerEdgeColor     	= 'none';
    sc.MarkerFaceColor    	= [dcol{2}];
    sc.MarkerFaceAlpha     	= alph;
    sc.MarkerFaceAlpha     	= alph;
end

ax5.XLabel.String           = 'Accuracy [norm]';
ax5.YLabel.String           = 'Lag [ms]';
ax5.FontSize                = lb_fs;
ax5.YLim                    = [400 1000];

lg                          = legend('','Faster','Slower');
lg.Location                 = 'north';
lg.FontSize                 = lg_fs;
lg.Box                      = 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AUROC solo - dyadic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

col                                 = jet(size(cr_dyad,1));
snr                                 = solo_perf(end).carr;

ax2                                 = axes('Position', [clm(1) height(2) dim]); hold on
pt                                  = patch([0 8 8 0], [.5 .5 1 1], [.75 1 .75], 'FaceAlpha',.1, 'EdgeColor','none');
pt                                  = patch([0 8 8 0], [0 0 .5 .5], [1 .75 .75], 'FaceAlpha',.1, 'EdgeColor','none');
tx                                  = text(1.25,.9, 'Dyad > Solo', 'FontSize', lb_fs, 'Color', 'k');
tx                                  = text(1.25,.1, 'Dyad < Solo', 'FontSize', lb_fs, 'Color', 'k');

ax3                                 = axes('Position', [clm(2) height(2) dim]); hold on
pt                                  = patch([0 8 8 0], [.5 .5 1 1], [.75 1 .75], 'FaceAlpha',.1, 'EdgeColor','none');
pt                                  = patch([0 8 8 0], [0 0 .5 .5], [1 .75 .75], 'FaceAlpha',.1, 'EdgeColor','none');

ax4                                 = axes('Position', [clm(3) height(2) dim]); hold on
pt                                  = patch([0 8 8 0], [.5 .5 1 1], [.75 1 .75], 'FaceAlpha',.1, 'EdgeColor','none');
pt                                  = patch([0 8 8 0], [0 0 .5 .5], [1 .75 .75], 'FaceAlpha',.1, 'EdgeColor','none');

for iD = 1:size(cr_dyad,1)
    
    if sum(contains(id_dyad{iD}, 'roh'))>0
        continue
    end
    
    % Extract performance data for each player in dyad
    for iP = 1:2
        clear d
        
        d.id                        = convertStringsToChars(id_dyad{iD}(iP));
        d.idx                       = find(strcmp(cellfun(@lower, sbj_lst, 'UniformOutput', false),d.id));
        d.sp                        = solo_perf(d.idx);
        d.sc                        = cr_solo{d.idx};
        d.dp                        = dyad_perf{iD,iP};
        d.dc                        = cr_dyad{iD,iP};
        
        for iCoh = 1:length(snr) 
            auc_acc(iD,iP,iCoh)   	= f_auroc(d.sp.acc{iCoh},d.dp.acc{iCoh});
            auc_ecc(iD,iP,iCoh)   	= f_auroc(d.sp.str{iCoh},d.dp.str{iCoh});
            auc_score(iD,iP,iCoh)  	= f_auroc(d.sp.trg_score{iCoh},d.dp.trg_score{iCoh});
            
            auc_cc(iD,iP,iCoh)    	= f_auroc(d.sc.cc(snr(iCoh) == d.sc.coh),d.dc.cc(snr(iCoh) == d.dc.coh));
            auc_xcp(iD,iP,iCoh)   	= f_auroc(d.sc.posPk(snr(iCoh) == d.sc.coh),d.dc.posPk(snr(iCoh) == d.dc.coh));
            auc_xc(iD,iP,iCoh)    	= f_auroc(d.sc.maxR(snr(iCoh) == d.sc.coh),d.dc.maxR(snr(iCoh) == d.dc.coh));
        end
    end
end

ax2                         = plotAUROC(ax2,auc_acc,{'Accuracy', '[AUROC]'},lb_fs,snr,alp,lw);
ax3                         = plotAUROC(ax3,auc_ecc,{'Eccentricity', '[AUROC]'},lb_fs,snr,alp,lw);
ax4                         = plotAUROC(ax4,auc_score,{'Score', '[AUROC]'},lb_fs,snr,alp,lw);

plotLines(ax2,lw)
plotLines(ax3,lw)
plotLines(ax4,lw)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AUROC summary - bar graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clm                                 = linspace(.2,.7,7);
dim                               	= [.075 .15];

% Individual player
for iCoh = 1:length(snr)
    ax                              = axes('Position', [clm(iCoh) height(3) dim]); hold on
    ax                              = plotBarGraphs(ax,auc_ecc, auc_acc, auc_score,iCoh,true,lb_fs);
    
    if iCoh == 1
        ax.YLabel.String            = 'Fraction [%]';
        ax.Title.String             = ['Coh = ' num2str(round(snr(iCoh),2))];
        
        lg                       	= legend('Better','Worse','Unchanged');
        lg.FontSize                 = lg_fs; 
        lg.Position(1:2)            = [.8 .205];
    else
        ax.YAxis.Visible            = 'off';
        ax.Title.String             = num2str(round(snr(iCoh),2));
    end
    
    ax.XAxis.Visible                = 'off';
    ax.Title.FontSize               = lb_fs;
end

%%% Display dyad-wise effects
for iCoh = 1:length(snr)
    ax                              = axes('Position', [clm(iCoh) height(4) dim]); hold on
    ax                              = plotBarGraphs(ax,auc_ecc, auc_acc, auc_score,iCoh,false,lb_fs);
  
    if iCoh == 1
        ax.YLabel.String            = 'Fraction [%]';
    else
        ax.XAxis.Visible            = 'off';
        ax.YAxis.Visible            = 'off';
    end
end

%%% Display leader-other effects

%%% Timing - Score correlation [Leader sorted]

%%% Look for oscillations in eccentricity signals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Annotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax0                         = axes('Position',[0 0 1 1],'Visible','off');
text(.03,height(3)+.075, {'Individual', 'Player'}, 'Parent', ax0, 'FontSize', 16, 'Color', 'k')
text(.03,height(4)+.075, {'Dyad'}, 'Parent', ax0, 'FontSize', 16, 'Color', 'k')

print(f, '/Users/fschneider/ownCloud/Documents/Publications/CPR_psychophysics/Figures/FIG4/FIG4', '-r400', '-dpng');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [solo_score, agnt_score, dyad_score] = getFinalScore(pth, sbj_lst)

cnt                                     = 0;

for iSubj = 1:length(sbj_lst)
    
    disp(['Score extraction subject: ' sbj_lst{iSubj}])
    data_pth                            = [pth sbj_lst{iSubj} '/summary/'];
    cd(data_pth)
    
    cnt                                 = cnt+1;
    mat_files                           = dir('*.mat');
    
    for iFile = 1:length(mat_files)
        clear tmp_tbl
        
        if contains(mat_files(iFile).name,'CPRsolo')
            tmp_tbl                     = load(mat_files(iFile).name);
            [trg_score{cnt,iFile}]      = getTargetScore(tmp_tbl.t);
            flag(cnt,iFile)             = 1;
            
        elseif contains(mat_files(iFile).name,'CPRagent')
            tmp_tbl                     = load(mat_files(iFile).name);
            
            if contains(mat_files(iFile).name,'agnt')
                continue
            end
            
            [trg_score{cnt,iFile}]      = getTargetScore(tmp_tbl.t);
            flag(cnt,iFile)             = 2;
        else
            [trg_score{cnt,iFile}]      = nan;
            flag(cnt,iFile)             = nan;
        end
    end
    
    dcnt = 0;
    
    for iDyad = 19:40
        cd(['/Volumes/T7_Shield/CPR_psychophysics/Dyad' num2str(iDyad) '/summary/'])
        
        dfiles                          = dir('*.mat');
        dcnt                            = dcnt+1;
        
        if isempty(dfiles)
            [dtrg_score{cnt,dcnt}]      = nan;
            dflag(cnt,dcnt)             = nan;
            continue
        end
        
        dnames                          = {dfiles.name};
        fsplit                          = cellfun(@(x) strsplit(x,'_'),dnames,'UniformOutput',false);
        subj                            = cellfun(@(x) x(2), fsplit);
        idx                             = cellfun(@(x) strcmp(x,lower(sbj_lst{iSubj})), subj);
        
        if sum(idx) == 0
            [dtrg_score{cnt,dcnt}]      = nan;
            dflag(cnt,dcnt)             = nan;
            continue
        else
            a                           = dnames(idx);
            b1                          = load(a{1});
            [dtrg_score{cnt,dcnt}]      = getTargetScore(b1.t);
            dflag(cnt,dcnt)             = 3;
            
            dcnt                    	= dcnt+1;
            b2                          = load(a{2});
            [dtrg_score{cnt,dcnt}]      = getTargetScore(b2.t);
            dflag(cnt,dcnt)             = 3;
        end
    end
    
    tmp_solo                                = trg_score(cnt,flag(cnt,:) == 1);
    tmp_solo(cellfun(@isempty, tmp_solo))   = [];
    
    tmp_agnt                                = trg_score(cnt,flag(cnt,:) == 2);
    tmp_agnt(cellfun(@isempty, tmp_agnt))   = [];
    
    tmp_dyad                                = dtrg_score(cnt,dflag(cnt,:) == 3);
    tmp_dyad(cellfun(@isempty, tmp_dyad))   = [];
    
    solo_score{cnt}                         = cellfun(@(x) x(end), tmp_solo);
    agnt_score{cnt}                         = cellfun(@(x) x(end), tmp_agnt);
    dyad_score{cnt}                         = cellfun(@(x) x(end), tmp_dyad);
end
end

function [trg_score] = getTargetScore(in)
c                           = 0;

% Target response params
for iState = 1:size(in.trg_ts,1)
    for iTrg = 1:length(in.trg_ts{iState})
        c                   = c +1;
        score_cum(c)        = in.trg_score{iState}(iTrg);
        score_coh(c)        = in.rdp_coh(iState);
        score_hi(c)         = in.trg_hit{iState}(iTrg);
    end
end
trg_score                   = score_cum(logical(score_hi));
end

function [tmp, tc] = extractTrials(trl_break, tbl, tmp, tc)
for iTrl = 2:length(trl_break) % modified function!
    tidx                    = trl_break(iTrl-1)+1: trl_break(iTrl);
    tt                      = tbl(tidx,:);
    
    % Initiate cell
    tc                      = tc+1;                             % Trial counter
    tmp.frme_ts{tc}         = [];
    tmp.rdp_dir{tc}         = [];
    tmp.rdp_coh{tc}     	= [];
    tmp.js_dir{tc}      	= [];
    tmp.js_str{tc}      	= [];
    tmp.refresh{tc}         = [];
    
    % Extract and add data
    for iState = 1:size(tt,1)
        tmp.frme_ts{tc}  	= [tmp.frme_ts{tc} tt.frme_ts{iState}];
        tmp.rdp_dir{tc}  	= [tmp.rdp_dir{tc} repmat(tt.rdp_dir(iState),1,length(tt.frme_ts{iState}))];
        tmp.rdp_coh{tc}  	= [tmp.rdp_coh{tc} repmat(tt.rdp_coh(iState),1,length(tt.frme_ts{iState}))];
        tmp.js_dir{tc}  	= [tmp.js_dir{tc} tt.js_dir{iState}];
        tmp.js_str{tc}  	= [tmp.js_str{tc} tt.js_str{iState}];
        tmp.refresh{tc}     = [tmp.refresh{tc} median(diff(tmp.frme_ts{tc}))];
    end
end
end

function out = response_readout(in, nSample)

c                           = 0;
snr                         = unique(in.rdp_coh);

% Target score
clear tscore score_cum score_hi score_coh trg_states
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

trg_states                 	= in.trg_hit(logical(in.trg_shown));
out.hir_pool                = sum(cellfun(@sum,trg_states)) / sum(cellfun(@numel,trg_states));

for iCoh = 1:length(snr)
    clear cIdx tIdx nhi ntrg
    
    cIdx = in.rdp_coh == snr(iCoh);
    tIdx = logical(in.trg_shown);
    
    tIdx(cellfun(@length,in.js_str) < 100) = false;
    cIdx(cellfun(@length,in.js_str) < 100) = false;
    
    % Hit rate
    nhi                     = sum(cellfun(@sum,in.trg_hit(cIdx & tIdx)));
    ntrg                    = sum(cellfun(@numel,in.trg_hit(cIdx & tIdx)));
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
    
    clear js_acc
    for iState = 1:length(rdp_dir)
        clear js_dev
        smpl_idx            = find(frmes{iState} < trg1_ts(iState),1,'last')-nSample : find(frmes{iState} < trg1_ts(iState),1,'last');
        js_dev              = rad2deg(circ_dist(deg2rad(js_dir{iState}(smpl_idx)),deg2rad(rdp_dir(iState))));  % Minimum RDP-Joystick difference
        js_acc(iState)      = nanmean(abs(1 - abs(js_dev) / 180));         	% Joystick accuracy
    end
    
    out.carr(iCoh)         	= snr(iCoh);
    out.macc(iCoh)        	= nanmedian(js_acc);
    out.acc{iCoh}        	= js_acc;
end
end

function ax = plotAUROC(ax,dat,str,lb_fs,snr,alp,lw)

axes(ax); hold on

% ln                          = line([0 8],[.5 .5], 'Color', 'k','LineStyle', '--', 'LineWidth', lw/2);

mat                         = [squeeze(dat(:,1,:)); squeeze(dat(:,2,:))];
mat(sum(mat,2) == 0,:)      = [];

for iL = 1:size(mat,1)
    pl(iL)               	= plot(mat(iL,:));
    pl(iL).Color          	= [.5 .5 .5 alp];
    pl(iL).LineWidth      	= lw/2;
end

pm                          = plot(nanmedian(mat),'LineWidth', lw, 'Color', [0 0 0]);

ax.XLim                  	= [1 length(snr)];
ax.YLim                  	= [0 1];
ax.FontSize               	= lb_fs;
ax.YLabel.String           	= str;
ax.XLabel.String           	= 'Coherence [%]';
ax.XTick                 	= 1:7;
ax.YTick                  	= [.2 .5 .8];
ax.XTickLabel            	= round(snr,2).*100;
ax.XTickLabelRotation      	= 0;
end

function plotLines(ax,lw)
axes(ax); hold on
n                           = .05;
crit                        = [.5-n .5+n];
lm                          = line([1 7],[.5 .5], 'Color', 'k', 'LineStyle', '--', 'LineWidth',lw);
ll                          = line([1 7],[crit(1) crit(1)], 'Color', 'k', 'LineStyle', '-.', 'LineWidth',lw/2);
lh                          = line([1 7],[crit(2) crit(2)], 'Color', 'k', 'LineStyle', '-.', 'LineWidth',lw/2);
end

function out = classifyEffects(mat,idx,flag)
cnt                         = 0;
n                           = .05;
crit                        = [.5-n .5+n];

if flag == true
    %%% Individual player
    for iD = 1:size(mat,1)
        for iP = 1:size(mat,2)
            cnt           	= cnt+1;
            better(cnt)   	= sum(mat(iD,iP,idx)>crit(2)) == length(mat(iD,iP,idx));
            worse(cnt)    	= sum(mat(iD,iP,idx)<crit(1)) == length(mat(iD,iP,idx));
            unclear(cnt)   	= better(cnt) == 0 && worse(cnt) == 0;
        end
    end 
else
    %%% Dyad-wise
    for iD = 1:size(mat,1)
        better(iD)       	= sum(mat(iD,1,idx)>crit(2)) == length(mat(iD,1,idx)) && sum(mat(iD,2,idx)>crit(2)) == length(mat(iD,2,idx));
        worse(iD)        	= sum(mat(iD,1,idx)<crit(1)) == length(mat(iD,1,idx)) && sum(mat(iD,2,idx)<crit(1)) == length(mat(iD,2,idx));
        unclear(iD)         = better(iD) == 0 && worse(iD) == 0;
    end
end

out                         = [sum(better) / numel(better)...
    sum(worse) / numel(worse)...
    sum(unclear) / numel(unclear)];
end

function ax = plotBarGraphs(ax,auc_ecc, auc_acc, auc_score,indx,plyer_flag,lb_fs)
bcol                                = [1 .5 .5; 0 0 0; .5 .5 .5];
a                                   = classifyEffects(auc_acc,indx,plyer_flag);
b                                   = classifyEffects(auc_ecc,indx,plyer_flag);
c                                   = classifyEffects(auc_score,indx,plyer_flag);
br                                  = bar([a;b;c],'stacked');

for i = 1:3
    br(i).FaceColor                 = bcol(i,:);
end

ax.XTick                            = 1:3;
ax.YTick                            = [0 .2 .4 .6 .8];
ax.XTickLabel                       = {'Accuracy','Eccentricity','Score'};
ax.YLim                             = [0 1];
ax.XLim                             = [0.25 3.75];
ax.FontSize                         = lb_fs;
ax.XAxis.FontSize                   = 11;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% UNUSED CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: RT vs XC %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% for iD = 1:size(cr_dyad,1)
%     %%% Timing leader
%     lgs(iD,:)             	= [median(cr_dyad{iD,1}.lag) median(cr_dyad{iD,2}.lag)];
%     [lgs(iD,:), plyer(iD,:)]= sort(lgs(iD,:));
%     
%     %%% Performance leader
%     if (nanmean(cell2mat(dyad_perf{iD,2}.trg_score)) - nanmean(cell2mat(dyad_perf{iD,1}.trg_score))) < 0
%         score_lead(iD,:)   	= [1 2];
%     elseif (nanmean(cell2mat(dyad_perf{iD,2}.trg_score)) - nanmean(cell2mat(dyad_perf{iD,1}.trg_score))) > 0
%         score_lead(iD,:)  	= [2 1];
%     end
%     
%     for iS = 1:2
%         clear rt
%         load([pth convertStringsToChars(id_dyad{iD}(plyer(iD,iS))) '/summary/' convertStringsToChars(id_dyad{iD}(plyer(iD,iS))) '_RT.mat'])
%         avg_rt(iD,iS)    	= median(rt.dat);
%         
%         try
%             sIdx                = cellfun(@(x) strcmp(convertStringsToChars(id_dyad{iD}(plyer(iD,iS))),x),cellfun(@lower,id_solo,'uni',false));
%             afc4_tresh(iD,iS)   = thresh_snr(sIdx);
%             avg_acc(iD,iS)      = nanmean(cell2mat(solo_perf(sIdx).acc));
%         catch
%             continue
%         end
%     end
% end
% 
% ax1                      	= axes('Position', [clm(2) height(1) dim]); hold on
% col                         = cbrewer('div', 'Spectral', size(cr_dyad,1), 'PCHIP');
% 
% sc                          = scatter(avg_rt(:,1) - avg_rt(:,2) , lgs(:,1) - lgs(:,2));
% sc.MarkerEdgeColor          = 'none';
% sc.MarkerFaceColor          = [0 0 0];
% 
% ax1.XLabel.String           = {'RT difference', '[Leader-Other]'};
% ax1.YLabel.String           = {'Lag difference', '[Leader-Other]'};
% % ax1.XLim                    = [-100 150];
% % ax1.YLim                    = [0 250];
% ax1.FontSize                = lb_fs;
% 
% ln                          = lsline;
% ln.LineWidth                = lw/2;
% ln.Color                    = [.75 .75 .75];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Dff Psychometric treshold vs Dff CPR accuracy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% avg_acc(11:end,:) = [];
% afc4_tresh(11:end,:) = [];
% 
% ax5                      	= axes('Position', [clm(3) height(1) dim]); hold on
% 
% sc                          = scatter(avg_acc(:,2) - avg_acc(:,1) , afc4_tresh(:,2) - afc4_tresh(:,1));
% sc.MarkerEdgeColor          = 'none';
% sc.MarkerFaceColor          = [0 0 0];
% 
% ax5.XLabel.String           = {'Accuracy difference', '[Leader-Other]'};
% ax5.YLabel.String           = {'4AFC thresh. difference', '[Leader-Other]'};
% % ax5.XLim                    = [-100 150];
% % ax5.YLim                    = [0 250];
% ax5.FontSize                = lb_fs;
% 
% ln                          = lsline;
% ln.LineWidth                = lw/2;
% ln.Color                    = [.75 .75 .75];
% Add relevant directories
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/WIP/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/Helper_functions/
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab
addpath /Users/fschneider/Documents/MATLAB/cbrewer/
addpath /Users/fschneider/Documents/MATLAB/palamedes1_10_9/Palamedes

close all
clear all

% Import subject summary spreadsheet
pth                         = '/Volumes/T7_Shield/CPR_psychophysics/';      % Local hard drive
x                           = readtable([pth 'Subjects_summary.xlsx']);     % Spreadsheet
sbj_lst                     = x.Abbreviation;                               % Subject ID list
sbj_lst(cellfun(@isempty,sbj_lst)) = [];

nSample                     = 30;                                           % Time window size [samples]
nLag                        = 150;                                          % Cross-correlation lag
s_cnt                       = 0;                                            % Subject counter
d_cnt                       = 0;                                            % Dyad counter

%% Extract solo data

% For all subjects
for iSubj = 1:length(sbj_lst)
    
    disp(['Processing subject: ' sbj_lst{iSubj}])
    
    data_pth                     	= [pth sbj_lst{iSubj} '/summary/'];   	% Data path
    solo_tbl                     	= [];                                   % Initialise
    dyad_tbl                     	= [];
    tmp_solo                        = [];
    tmp_dyad                        = [];
    sc_cnt                       	= 0;                                    % Cycle counter
    dc_cnt                       	= 0;
    
    if isdir(data_pth)
        cd(data_pth)
        mat_files               	= dir('*.mat');
        
        % For all files in directory
        for iFile = 1:length(mat_files)
            
            if contains(mat_files(iFile).name,'CPRsolo')
                % Load pre-processed data table
                tmp_tbl             = load(mat_files(iFile).name);
                
                % Exclude sessions with different coherence level
                if sum(unique(tmp_tbl.t.rdp_coh) < .1) > 2
                    continue
                end
                
                % Concatenate session summary tables
                solo_tbl           	= [solo_tbl; tmp_tbl.t];                % Organise solo data
            end
        end
        
        % For all dyads
        dyad_cnt                    = 0;
        for iDyad = 19:61
            cd([pth ['Dyad' num2str(iDyad)] '/summary/'])
            mat_files              	= dir('*.mat');
            
            if isempty(mat_files)
                continue
            end
            dyad_cnt                = dyad_cnt+1;
            tmp1                    = split(mat_files(1).name,'_');
            tmp3                    = split(mat_files(3).name,'_');
            id_dyad(dyad_cnt,:)    	= [{tmp1{2}}, {tmp3{2}}];
            n_dyad(dyad_cnt,:)      = iDyad;
            
            % Check if subject part of this dyad
            for iFile = 1:size(mat_files,1)
                is_player(iFile) = contains(mat_files(iFile).name,lower(sbj_lst{iSubj}));
                
                if is_player(iFile)
                    tmp_tbl      	= load(mat_files(iFile).name);
                    dyad_tbl      	= [dyad_tbl; tmp_tbl.t];                % Organise dyadic data
                end
            end
        end
        
        if ~isempty(solo_tbl)
            % Performance [time window analysis]
            solo_perf{iSubj}      	= response_readout(solo_tbl, nSample);

            % Extract stimulus cycles for correlation analysis
            cyc_boundary        	= [0; find(diff(solo_tbl.cyc_no) ~= 0)];
            [tmp_solo, sc_cnt]   	= extractCycles(cyc_boundary, solo_tbl, tmp_solo, sc_cnt);
            
            % Correlation analysis
            [solo_cr{iSubj},ps]     = CPR_correlation_analysis_WIP(tmp_solo, nLag, false);
            
            solo_perf{iSubj}.id     = lower(sbj_lst{iSubj});
            solo_cr{iSubj}.id       = lower(sbj_lst{iSubj});
            
            
        end
        
        if ~isempty(dyad_tbl)
            % Performance [time window analysis]
            dyad_perf{iSubj}      	= response_readout(dyad_tbl, nSample);
            
            % Extract stimulus cycles for correlation analysis
            cyc_boundary        	= [0; find(diff(dyad_tbl.cyc_no) ~= 0)];
            [tmp_dyad, dc_cnt]   	= extractCycles(cyc_boundary, dyad_tbl, tmp_dyad, dc_cnt);
            
            % Correlation analysis
            [dyad_cr{iSubj},ps]     = CPR_correlation_analysis_WIP(tmp_dyad, nLag, false);
            
            dyad_perf{iSubj}.id     = lower(sbj_lst{iSubj});
            dyad_cr{iSubj}.id       = lower(sbj_lst{iSubj});
        end
    end
end

%% Convert relevant data to matrix for simplicity

cnt                         = 0;
for iSubj = 1:length(sbj_lst)
    if isempty(solo_perf{iSubj}.hir)
        continue
    end
    
    cnt                   	= cnt + 1;
    solo_hir(cnt,:)        	= solo_perf{iSubj}.hir;          % Hit rate
    solo_trg_score(cnt,:)  	= solo_perf{iSubj}.trg_mscore;   % Avg target scores
    solo_macc(cnt,:)       	= solo_perf{iSubj}.macc_trg;     % Avg accuracy
    solo_mecc(cnt,:)       	= solo_perf{iSubj}.mecc;         % Avg eccentricity
    
    try
        dyad_hir(cnt,:)       	= dyad_perf{iSubj}.hir;          % Hit rate
        dyad_trg_score(cnt,:)	= dyad_perf{iSubj}.trg_mscore;   % Avg target scores
        dyad_macc(cnt,:)     	= dyad_perf{iSubj}.macc_trg;     % Avg accuracy
        dyad_mecc(cnt,:)    	= dyad_perf{iSubj}.mecc;         % Avg eccentricity
    catch
        dyad_hir(cnt,:)       	= nan;
        dyad_trg_score(cnt,:)	= nan;
        dyad_macc(cnt,:)     	= nan;
        dyad_mecc(cnt,:)    	= nan;
    end
end

%% PLOT

% f                           = figure('units','normalized','position',[0 0 .5 1]);
f                           = figure('units','centimeters','position',[0 0 21 29.7]);
height                      = [linspace(.555, .056,4) .1];
colmn                       = linspace(.075, .85,4);
lb_fs                       = 14;
lg_fs                       = 10;
lw                          = 3;
frme_ms                     = 1000/120;
alp                         = .4;
dim                         = [.18 .15];
col_dat                     = [0 0 0];
col_ci                      = [.3 0 0];
snr                         = solo_perf{end}.carr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Score comparison %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sbjs = cellfun(@lower, sbj_lst, 'UniformOutput', false);

pl_dim                      = dim*1.33;
row                         = .79;
cmap_coh                    = cool(size(snr,2));

ax20                    	= axes('Position', [colmn(1) row pl_dim]); hold on

for iSubj = 1:length(dyad_perf)
    if isempty(dyad_perf{iSubj})
        continue
    end
    sidx = cellfun(@(x) strcmp(x,dyad_perf{iSubj}.id),sbjs);
    for iCoh = 1:length(snr)
        sc(iCoh)                    = scatter(solo_perf{sidx}.trg_mscore(iCoh),dyad_perf{iSubj}.trg_mscore(iCoh));
        sc(iCoh) .MarkerFaceColor   = cmap_coh(iCoh,:);
        sc(iCoh) .MarkerFaceAlpha   = alp;
        sc(iCoh) .MarkerEdgeColor   = 'none';
        lg_str{iCoh}               	= num2str(round(snr(iCoh),2)*100);
    end
end

ln                          = line([0 1],[0 1]);
ln.LineStyle                = ':';
ln.Color                    = [0 0 0];

ax20.FontSize               = lb_fs;
ax20.YLabel.String          = 'Dyad score [a.u.]';
ax20.XLabel.String          = 'Solo score [a.u.]';
ax20.XTickLabelRotation     = 0;

lg0                         = legend(sc,lg_str,'Location','northwest','NumColumns', 1);
lg0.Box                     = 'off';
lg0.Position(1)             = .4;

% Implement gramm solution later!
% clear g
% g=gramm('x',xsolo,'y',ydyad,'color',[0 0 0]);
% g.facet_grid([],b);
% g.geom_point();
% g.stat_cornerhist('aspect',0.5);
% g.geom_abline();
% g.set_title('Visualize x-y with stat_cornerhist()');
% figure('Position',[100 100 800 600]);
% g.draw();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Cumulative score comparison %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cnt = 0;
for iSubj = 1:length(dyad_perf)
    if isempty(dyad_perf{iSubj})
        continue
    end
    
    cnt                     = cnt+1;
    sidx                    = cellfun(@(x) strcmp(x,dyad_perf{iSubj}.id),sbjs);        
    mscore_solo(cnt)    	= mean(solo_perf{sidx}.score_norm);
    mscore_dyadic(cnt)   	= mean(dyad_perf{iSubj}.score_norm);
end

ax21                    	= axes('Position', [colmn(2)+.05 row pl_dim]); hold on

ln                          = line([0 1],[0 1]);
ln.LineStyle                = ':';
ln.Color                    = [0 0 0];

sc                          = scatter(mscore_solo,mscore_dyadic);
sc.MarkerFaceColor          = [.3 .3 .3];
sc.MarkerFaceAlpha          = 1;
sc.MarkerEdgeColor          = 'none';

ax21.FontSize               = lb_fs;
ax21.YLim                   = [0 .4];
ax21.XLim                   = [0 .4];
ax21.YTick                  = [0 .2 .4];
ax21.XTick                  = [0 .2 .4];
ax21.YLabel.String          = 'Dyad score [a.u.]';
ax21.XLabel.String          = 'Solo score [a.u.]';
ax21.XTickLabelRotation     = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: SOLO - Hit rate raw %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax1                         = axes('Position', [colmn(1) height(1) dim]); hold on
[ax1,pl]                    = plotData(ax1,solo_hir,true,lw,alp,col_dat,col_ci);
ax1.YLim                    = [10 70];
ax1.XLim                    = [1 size(solo_hir,2)];
ax1.XLabel.String           = 'Coherence [%]';
ax1.YLabel.String           = 'Hit rate [%]';
ax1.YLabel.Position(1)      = -.4;
ax1.XTick                   = 1:length(snr);
ax1.XTickLabel              = round(snr,2)*100;
ax1.FontSize                = lb_fs;
ax1.XTickLabelRotation      = 0;
ax1.XAxis.Visible           = 'off';
% ax1.XGrid                   = 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: SOLO - Lag %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iSubj = 1:size(solo_cr,2)
    solo_mlag_id{iSubj} = solo_cr{iSubj}.id;

    for iCoh = 1:length(snr)
        clear cIdx
        cIdx             	= solo_cr{iSubj}.coh == snr(iCoh);
        solo_mlag(iSubj,iCoh) = mean(solo_cr{iSubj}.lag(cIdx));
    end
end

ax2                         = axes('Position', [colmn(1) height(4) dim]); hold on
[ax2,pl]                    = plotData(ax2,solo_mlag ./ 1e3,false,lw,alp,col_dat,col_ci);
ax2.XLim                    = [1 size(solo_mlag,2)];
ax2.YLim                    = [.3 1];
ax2.XLabel.String           = 'Coherence [%]';
ax2.YLabel.String           = 'Lag [s]';
ax2.YLabel.Position(1)      = -.4;
ax2.XTick                   = 1:length(snr);
ax2.XTickLabel              = round(snr,2)*100;
ax2.FontSize                = lb_fs;
ax2.XTickLabelRotation      = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: SOLO - Avg accuracy raw %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax3                         = axes('Position', [colmn(1) height(2) dim]); hold on
[ax3,pl]                    = plotData(ax3,solo_macc,true,lw,alp,col_dat,col_ci);
ax3.XLim                    = [1 size(solo_macc,2)];
ax3.YLim                    = [30 100];
ax3.XLabel.String           = 'Coherence [%]';
ax3.YLabel.String           = 'Accuracy [%]';
ax3.YLabel.Position(1)      = -.4;
ax3.XTick                   = 1:length(snr);
ax3.XTickLabel              = round(snr,2)*100;
ax3.FontSize                = lb_fs;
ax3.XTickLabelRotation      = 0;
ax3.XAxis.Visible           = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: SOLO - Avg eccentricity raw %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax4                         = axes('Position', [colmn(1) height(3) dim]); hold on
[ax4,pl]                    = plotData(ax4,solo_mecc,true,lw,alp,col_dat,col_ci);
ax4.XLim                    = [1 size(solo_mecc,2)];
ax4.YLim                    = [20 100];
ax4.XLabel.String           = 'Coherence [%]';
ax4.YLabel.String           = 'Eccentricity [%]';
ax4.YLabel.Position(1)      = -.4;
ax4.XTick                   = 1:length(snr);
ax4.XTickLabel              = round(snr,2)*100;
ax4.FontSize                = lb_fs;
ax4.XTickLabelRotation      = 0;
ax4.XAxis.Visible           = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: DYADIC - Hit rate raw %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax5                         = axes('Position', [colmn(2) height(1) dim]); hold on
[ax5]                       = plotScatter(ax5, solo_hir, dyad_hir, snr, lb_fs, true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Dyadic - Lag %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub_cnt = 0;
for iSubj = 1:size(dyad_cr,2)
    if isempty(dyad_cr{iSubj})
        continue
    end
    
    sub_cnt                     = sub_cnt+1;
    dyad_mlag_id{sub_cnt}     	= dyad_cr{iSubj}.id;
        
    for iCoh = 1:length(snr)
        clear cIdx
        cIdx                	= dyad_cr{iSubj}.coh == snr(iCoh);
        dyad_mlag(sub_cnt,iCoh)	= mean(dyad_cr{iSubj}.lag(cIdx));
    end
end

for iSubj = 1:size(dyad_mlag_id,2)
    idx                         = cellfun(@(x) strcmp(x,dyad_mlag_id{iSubj}),solo_mlag_id);
    solo_mlag_sorted(iSubj,:)   = solo_mlag(idx,:);
end

ax6                         = axes('Position', [colmn(2) height(4) dim]); hold on
ln                          = line([0 1000],[0 1000]);
ln.LineStyle                = ':';
ln.Color                    = [0 0 0];
[ax6]                       = plotScatter(ax6, solo_mlag_sorted ./ 1e3, dyad_mlag ./ 1e3, snr, lb_fs, false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: DYADIC - Avg accuracy raw %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax7                         = axes('Position', [colmn(2) height(2) dim]); hold on
[ax7]                       = plotScatter(ax7, solo_macc, dyad_macc, snr, lb_fs, true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: DYADIC - Avg eccentricity raw %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax8                         = axes('Position', [colmn(2) height(3) dim]); hold on
[ax8]                       = plotScatter(ax8, solo_mecc, dyad_mecc, snr, lb_fs, true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBPLOT: Hit rate comparison %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cnt = 0;
for iSubj = 1:size(dyad_perf,2)
    if ~isempty(dyad_perf{iSubj})
        idx                         = cellfun(@(x) strcmp(x,dyad_perf{iSubj}.id),lower(sbj_lst));        
        cnt                         = cnt+1;
        df_hir(cnt,:)               = solo_perf{idx}.hir - dyad_perf{iSubj}.hir;
        df_hir_pool(cnt,:)          = solo_perf{idx}.hir_pool - dyad_perf{iSubj}.hir_pool;
    end
end

ax13                      	= axes('Position', [colmn(3) height(1) dim]); hold on
% pt                        = patch([0 8 8 0], [0 0 25 25], [.7 .7 .7], 'FaceAlpha',.2, 'EdgeColor','none');
% pt                       	= patch([0 8 8 0], [-25 -25 0 0], [.3 .3 .3], 'FaceAlpha',.2, 'EdgeColor','none');
tx                        	= text(1.25,18, 'Dyad > Solo', 'FontSize', lb_fs, 'Color', 'k');
tx                       	= text(1.25,-18, 'Dyad < Solo', 'FontSize', lb_fs, 'Color', 'k');

[ax13,pl]                   = plotData(ax13,df_hir,true,lw,alp,col_dat,col_ci);
lm                          = line([1 7],[0 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth',lw);
ax13.XLim                   = [1 size(df_hir,2)];
ax13.FontSize               = lb_fs;
ax13.YLabel.String          = 'Difference';
ax13.XAxis.Visible          = 'off';
ax13.YLim                   = [-25 25];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUROC solo - dyadic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:size(dyad_perf,2)
    if isempty(dyad_perf{j})
        ply_id{j} = nan;
    else
        ply_id{j} = dyad_perf{j}.id;
    end
end

cnt = 0;

for iSubj = 1:size(dyad_perf,2)
    
    d.sp                  	= solo_perf{iSubj};                             % Performance data
    d.sc                  	= solo_cr{iSubj};                               % Correlation data
    
    sIdx                    = cellfun(@(x) strcmp(x,solo_perf{iSubj}.id),ply_id);
    
    if sum(sIdx) == 0
        continue
    end
    
    d.dp                  	= dyad_perf{sIdx};
    d.dc                  	= dyad_cr{sIdx};
         
    if isempty(d.dc)
        continue
    end
    
    cnt = cnt+1;
    for iCoh = 1:length(snr)
        auc_acc(cnt,iCoh)       = f_auroc(d.sp.acc_trg{iCoh},d.dp.acc_trg{iCoh});
%         auc_acc(iScntubj,iCoh)  	= f_auroc(d.sp.acc_state{iCoh},d.dp.acc_state{iCoh});
        auc_ecc(cnt,iCoh)   	= f_auroc(d.sp.ecc{iCoh},d.dp.ecc{iCoh});
        auc_score(cnt,iCoh)  	= f_auroc(d.sp.trg_score{iCoh},d.dp.trg_score{iCoh});
        
        auc_cc(cnt,iCoh)    	= f_auroc(d.sc.cc(snr(iCoh) == d.sc.coh),d.dc.cc(snr(iCoh) == d.dc.coh));
        auc_xcp(cnt,iCoh)   	= f_auroc(d.sc.posPk(snr(iCoh) == d.sc.coh),d.dc.posPk(snr(iCoh) == d.dc.coh));
        auc_xc(cnt,iCoh)    	= f_auroc(d.sc.maxR(snr(iCoh) == d.sc.coh),d.dc.maxR(snr(iCoh) == d.dc.coh));
        auc_lag(cnt,iCoh)    	= f_auroc(d.sc.lag(snr(iCoh) == d.sc.coh),d.dc.lag(snr(iCoh) == d.dc.coh));
                
        p_acc(cnt,iCoh)      	= ranksum(d.sp.acc_trg{iCoh},d.dp.acc_trg{iCoh});
%         p_acc(cnt,iCoh)      	= ranksum(d.sp.acc_state{iCoh},d.dp.acc_state{iCoh});
        p_ecc(cnt,iCoh)      	= ranksum(d.sp.ecc{iCoh},d.dp.ecc{iCoh});
        p_score(cnt,iCoh)       = ranksum(d.sp.trg_score{iCoh},d.dp.trg_score{iCoh});
        p_cc(cnt,iCoh)      	= ranksum(d.sc.cc(snr(iCoh) == d.sc.coh),d.dc.cc(snr(iCoh) == d.dc.coh));
        p_xcp(cnt,iCoh)      	= ranksum(d.sc.posPk(snr(iCoh) == d.sc.coh),d.dc.posPk(snr(iCoh) == d.dc.coh));
        p_xc(cnt,iCoh)       	= ranksum(d.sc.maxR(snr(iCoh) == d.sc.coh),d.dc.maxR(snr(iCoh) == d.dc.coh));
        p_lag(cnt,iCoh)      	= ranksum(d.sc.lag(snr(iCoh) == d.sc.coh),d.dc.lag(snr(iCoh) == d.dc.coh));
    end
    
    p_ecc_pooled(cnt)           = ranksum(cell2mat(d.sp.ecc'),cell2mat(d.dp.ecc'));
    auc_ecc_pooled(cnt)         = f_auroc(cell2mat(d.sp.ecc'),cell2mat(d.dp.ecc'));
    
    p_acc_pooled(cnt)           = ranksum(cell2mat(d.sp.acc_trg),cell2mat(d.dp.acc_trg));
    auc_acc_pooled(cnt)         = f_auroc(cell2mat(d.sp.acc_trg),cell2mat(d.dp.acc_trg));
end

ax9                       	= axes('Position', [colmn(3) height(4) dim]); hold on

ax10                     	= axes('Position', [colmn(3) height(2) dim]); hold on
tx                        	= text(1.25,.8, 'Dyad > Solo', 'FontSize', lb_fs, 'Color', 'k');
tx                       	= text(1.25,.2, 'Dyad < Solo', 'FontSize', lb_fs, 'Color', 'k');

ax11                       	= axes('Position', [colmn(3) height(3) dim]); hold on
% pt                         	= patch([0 8 8 0], [.5 .5 1 1], [.75 .75 .75], 'FaceAlpha',.1, 'EdgeColor','none');
% pt                         	= patch([0 8 8 0], [0 0 .5 .5], [.5 .5 .5], 'FaceAlpha',.1, 'EdgeColor','none');

ax9                     	= plotAUROC(ax9,auc_lag,'AUROC',lb_fs,snr,alp,lw,col_dat,col_ci);
ax9.XAxis.Visible         	= 'on';
ax10                     	= plotAUROC(ax10,auc_acc,'AUROC',lb_fs,snr,alp,lw,col_dat,col_ci);
ax11                    	= plotAUROC(ax11,auc_ecc,'AUROC',lb_fs,snr,alp,lw,col_dat,col_ci);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUROC Significant subjects histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nBin                        = 15;
ax12                       	= axes('Position', [colmn(3)+.08 row pl_dim]); hold on
hs                          = histogram(auc_ecc_pooled(p_ecc_pooled < .05/length(p_ecc_pooled)),nBin);
hs.FaceColor                = [.3 .3 .3];
hs.FaceAlpha                = 1;
ax12.FontSize               = lb_fs;
ax12.YLim                   = [0 5];
ax12.YLabel.String          = '# Subjects';
ax12.XLabel.String          = 'AUROC';
ln                          = line([.5 .5],[0 5]);
ln.LineWidth                = 1.5;
ln.LineStyle                = ':';
ln.Color                    = [0 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUROC Significance testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sig_boundary                = .05 / (size(p_ecc,1) * size(p_ecc,2)); % Bonferroni correction

ax15                       	= axes('Position', [colmn(4) height(2) .1 dim(2)]); hold on
ax15                      	= plotBar(ax15, auc_acc, p_acc, sig_boundary, lb_fs, snr);

ax16                       	= axes('Position', [colmn(4) height(3) .1 dim(2)]); hold on
ax16                      	= plotBar(ax16, auc_ecc, p_ecc, sig_boundary, lb_fs, snr);

ax17                       	= axes('Position', [colmn(4) height(4) .1 dim(2)]); hold on
ax17                      	= plotBar(ax17, auc_lag, p_lag, sig_boundary, lb_fs, snr);
ax17.XAxis.Visible         	= 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Annotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dest_dir                    = '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG2/';
ax00                         = axes('Position',[0 0 1 1],'Visible','off');

% vofs = .075;
% hofs = .16;
% text(colmn(1)-vofs,.98, 'A', 'Parent', ax00, 'FontSize', 22, 'Color', 'k')
% text(colmn(3)-.125,.98, 'B', 'Parent', ax00, 'FontSize', 22, 'Color', 'k')
% text(colmn(1)-vofs,height(1)+hofs, 'C', 'Parent', ax00, 'FontSize', 22, 'Color', 'k')
% text(colmn(1)-vofs,height(2)+hofs, 'D', 'Parent', ax00, 'FontSize', 22, 'Color', 'k')
% text(colmn(1)-vofs,height(3)+hofs, 'E', 'Parent', ax00, 'FontSize', 22, 'Color', 'k')
% text(colmn(1)-vofs,height(4)+hofs, 'F', 'Parent', ax00, 'FontSize', 22, 'Color', 'k')

ofs = 0;
text(colmn(1)+ofs,.71, 'Solo condition', 'Parent', ax00, 'FontSize', 16, 'Color', 'k')
text(colmn(2)+ofs,.71, 'Solo vs Dyadic', 'Parent', ax00, 'FontSize', 16, 'Color', 'k')
text(colmn(3)+ofs,.71, 'Effect size', 'Parent', ax00, 'FontSize', 16, 'Color', 'k')
text(colmn(4)+ofs,.55, 'Stats', 'Parent', ax00, 'FontSize', 16, 'Color', 'k')

print(f, [dest_dir '/FIG2'], '-r500', '-dpng');
print(f, [dest_dir '/FIG2'], '-r500', '-dsvg', '-painters');
print(f, [dest_dir '/FIG2'], '-r500', '-depsc2', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Reported stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Eccentricity [coherence pooled, within subject]
sign_bool      	= p_ecc_pooled < ( .05 / length(p_ecc_pooled));
perc_sign       = sum(sign_bool) / length(sign_bool);

% Accuracy [coherence pooled, within subject]
sign_bool      	= p_acc_pooled < ( .05 / length(p_acc_pooled));
perc_sign       = sum(sign_bool) / length(sign_bool);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tmp, cc] = extractCycles(cyc_break, tbl, tmp, cc)
for iCyc = 2:length(cyc_break)
    tidx                    = cyc_break(iCyc-1)+1 : cyc_break(iCyc);
    tt                      = tbl(tidx,:);
    
    % Initiate cell
    cc                      = cc+1;     % Cycle counter
    tmp.frme_ts{cc}         = [];
    tmp.rdp_dir{cc}         = [];
    tmp.rdp_coh{cc}     	= [];
    tmp.js_dir{cc}      	= [];
    tmp.js_ecc{cc}      	= [];
    tmp.refresh{cc}         = [];
    
    % Extract and add data
    for iState = 1:size(tt,1)
        tmp.frme_ts{cc}  	= [tmp.frme_ts{cc} tt.frme_ts{iState}];
        tmp.rdp_dir{cc}  	= [tmp.rdp_dir{cc} repmat(tt.rdp_dir(iState),1,length(tt.frme_ts{iState}))];
        tmp.rdp_coh{cc}  	= [tmp.rdp_coh{cc} repmat(tt.rdp_coh(iState),1,length(tt.frme_ts{iState}))];
        tmp.js_dir{cc}  	= [tmp.js_dir{cc} tt.js_dir{iState}];
        tmp.js_ecc{cc}  	= [tmp.js_ecc{cc} tt.js_ecc{iState}];
        tmp.refresh{cc}     = [tmp.refresh{cc} median(diff(tmp.frme_ts{cc}))];
    end
end
end

function out = response_readout(in, nSample)

c                           = 0;
snr                         = unique(in.rdp_coh);

% Target score
clear tscore score_cum score_hi score_coh score_dte score_exp trg_states
for iState = 1:size(in.trg_ts,1)
    for iTrg = 1:length(in.trg_ts{iState})
        try
            c                   = c +1;
            score_hi(c)         = in.trg_hit{iState}(iTrg);
            score(c)            = double(in.trg_score{iState}(iTrg));
            score_coh(c)        = in.rdp_coh(iState);
            score_dte(c)        = in.date(iState);
            score_exp(c)        = in.exp(iState);
        catch
            warning(['Skipped state/target: ' num2str(iState) '/' num2str(iTrg)])
        end
    end
end

trg_states                 	= in.trg_hit(logical(in.trg_shown));
out.hir_pool                = sum(cellfun(@sum,trg_states)) / sum(cellfun(@numel,trg_states));

dte = unique(score_dte);
dte(ismissing(dte)) = [];

for iDate = 1:length(dte)
    clear dte_idx score_cum dte_score
    dte_idx                 = score_dte == dte(iDate);
    out.score_final_exp(iDate) = unique(score_exp(dte_idx));
    dte_score               = score(dte_idx);
    score_cum               = cumsum(dte_score(~isnan(dte_score)));
    out.score_final(iDate)  = score_cum(end);
    out.score_norm(iDate)   = score_cum(end) ./ length(dte_score(~isnan(dte_score)));
end

for iCoh = 1:length(snr)
    clear cIdx tIdx nhi ntrg
    
    cIdx = in.rdp_coh == snr(iCoh);
    tIdx = logical(in.trg_shown);
    
    tIdx(cellfun(@length,in.js_ecc) < 100) = false;
    cIdx(cellfun(@length,in.js_ecc) < 100) = false;
    
    % Hit rate
    nhi                     = sum(cellfun(@sum,in.trg_hit(cIdx & tIdx)));
    ntrg                    = sum(cellfun(@numel,in.trg_hit(cIdx & tIdx)));
    out.hir(iCoh)           = nhi / ntrg;
    
    % Target score [hits only]
    out.trg_mscore(iCoh)	= nanmean(score(score_coh  == snr(iCoh) & score_hi == true));
    out.trg_score{iCoh}  	= score(score_coh  == snr(iCoh) & score_hi == true);
    
    % Joystick displacement
    out.mecc(iCoh)         	= nanmedian(cellfun(@(x) nanmedian(x(end-nSample:end)), in.js_ecc(cIdx)));
    out.ecc{iCoh}           = cellfun(@(x) nanmedian(x(end-nSample:end)), in.js_ecc(cIdx));
    
    % Joystick accuracy [state-wise]
    for iState = 1:length(in.rdp_dir)
        % At least 100 frames
        if length(in.js_dir{iState}) < 100
            continue
        end
        js_dev              = rad2deg(circ_dist(deg2rad(in.js_dir{iState}(end-nSample:end)),deg2rad(in.rdp_dir(iState))));  % Minimum RDP-Joystick difference
        js_acc(iState)      = nanmean(abs(1 - abs(js_dev) / 180));         	% Joystick accuracy
    end
    out.macc_state(iCoh)  	= nanmedian(js_acc);
    out.acc_state{iCoh}    	= js_acc;
    
    % Joystick accuracy [before first target]
    t1_ts                   = cellfun(@(x) x(1), in.trg_ts);
    f1_ts                   = cellfun(@(x) x(1), in.frme_ts);
    trgIdx                  = (t1_ts-f1_ts) >= 1e6;
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
    
    out.macc_trg(iCoh)      = nanmedian(js_acc);
    out.acc_trg{iCoh}     	= js_acc;
    out.carr(iCoh)         	= snr(iCoh);
end
end

function [ax,pl] = plotData(ax,dat,scale_flag,lw,alp,col_dat,col_ci)

% Remove NaNs and zero rows
dat(sum(isnan(dat),2)>0,:)  = [];
dat(sum(dat==0,2)>0,:)      = [];

if scale_flag
    dat                     = dat.*100;
end

% Boostrap confidence intervals
nRep                        = 1000;
[CI,~]                      = bootci(nRep,{@mean,dat},'Alpha',0.05);

% Prepare filled area
vec                         = 1:length(CI);
x_spacing                   = [vec fliplr(vec)];
ci                          = [CI(1,:) fliplr(CI(2,:))];

% Plot subject-wise data
hold on
for iL = 1:size(dat,1)
    pl(iL)                  = plot(dat(iL,:), 'LineWidth', lw/2, 'Color', [.5 .5 .5 alp]);
end

% Overlay confidence intervals
fl                          = fill(x_spacing,ci,col_ci,'EdgeColor','none', 'FaceAlpha', alp);

% Plot mean curve
pl                          = plot(mean(dat), 'LineWidth', lw, 'Color', col_dat);
end

function [ax] = plotScatter(ax, solo_dat, dyad_dat, snr, lb_fs, ax_flag)

hold on
coh_col                     = cool(length(snr));

for iCoh = 1:size(solo_dat,2)
    sc                      = scatter(solo_dat(:,iCoh),dyad_dat(:,iCoh));
    sc.MarkerFaceColor      = coh_col(iCoh,:);
    sc.MarkerEdgeColor      = 'none';
    sc.MarkerFaceAlpha      = .4;
end

ax.FontSize                 = lb_fs;
ax.XLabel.String            = 'Solo';
ax.YLabel.String            = 'Dyad';
ax.XLim                 	= [0 1];
ax.YLim                     = [0 1];
ax.XTick                    = [.1 .5 .9];
ax.YTick                    = [.1 .5 .9];

if ax_flag 
    ax.XAxis.Visible        = 'off';
end

ln                          = line([0 1],[0 1]);
ln.LineWidth                = 1.5;
ln.LineStyle                = ':';
ln.Color                    = [0 0 0];

end

function ax = plotAUROC(ax,dat,str,lb_fs,snr,alp,lw,col_dat,col_ci)

axes(ax); hold on

% Remove zero rows
dat(sum(dat == 0,2) == size(dat,2),:) = [];

% Plot subject-wise data
for iL = 1:size(dat,1)
    pl(iL)               	= plot(dat(iL,:));
    pl(iL).Color          	= [.5 .5 .5 alp];
    pl(iL).LineWidth      	= lw/2;
end

% Boostrap confidence intervals
nRep                        = 1000;
[CI,~]                      = bootci(nRep,{@mean,dat},'Alpha',0.05);

% Prepare filled area
vec                         = 1:length(CI);
x_spacing                   = [vec fliplr(vec)];
ci                          = [CI(1,:) fliplr(CI(2,:))];

% Overlay confidence intervals
fl                          = fill(x_spacing,ci,col_ci,'EdgeColor','none', 'FaceAlpha', alp);

% Plot mean curve
pm                          = plot(mean(dat),'LineWidth', lw, 'Color', col_dat);

ax.XLim                  	= [1 length(snr)];
ax.YLim                  	= [.1 .9];
ax.FontSize               	= lb_fs;
ax.YLabel.String           	= str;
ax.XLabel.String           	= 'Coherence [%]';
ax.XTick                 	= 1:7;
ax.YTick                  	= [.25 .5 .75];
ax.XTickLabel            	= round(snr,2).*100;
ax.XTickLabelRotation      	= 0;
ax.XAxis.Visible            = 'off';

% n                           = .05;
% crit                        = [.5-n .5+n];
lm                          = line([1 7],[.5 .5], 'Color', 'k', 'LineStyle', '--', 'LineWidth',lw);
% ll                          = line([1 7],[crit(1) crit(1)], 'Color', 'k', 'LineStyle', '-.', 'LineWidth',lw/2);
% lh                          = line([1 7],[crit(2) crit(2)], 'Color', 'k', 'LineStyle', '-.', 'LineWidth',lw/2);

end

function ax = plotBar(ax,dat,p_mat,p_val,lb_fs,snr)

axes(ax); hold on

pos = (dat > .5) & (p_mat < p_val);
neg = (dat < .5) & (p_mat < p_val);
all = size(dat,1);

bp                          = bar( (sum(pos)/all) .*100);
bp.FaceColor                = [.4 .4 .4];

bn                          = bar( -(sum(neg)/all) .*100);
bn.FaceColor                = [.1 .1 .1];

ax.YLim                     = [-60 60];
ax.YLabel.String            = '[%]';
ax.FontSize                 = lb_fs;
ax.XTick                 	= [1 4 7];
ax.XTickLabel            	= round(snr([1 4 7]),2).*100;
ax.XTickLabelRotation      	= 0;
ax.XAxis.Visible         	= 'off';
ax.YLabel.String           	= '[%]';
ax.XLabel.String           	= 'Coherence [%]';
end

%% STUFF NO LONGER NEEDED 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% SUBPLOT: Example subject joystick response %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% row                         = .79;
% clm                         = .15;
% pl_dim                      = dim*1.2;
% 
% cmap                        = [1 1 1; flipud(gray(256))];
% steps                       = .05;
% bins                        = 0:steps:1;
% c                           = 0;
% tbl                         = load([pth 'AnM/summary/20220629_anm_CPRsolo_block2_tbl.mat']);
% t                           = tbl.t;
% 
% ax0h                        = axes('Position', [clm-.082 row pl_dim(1)/5 pl_dim(2)]); hold on
% ax0v                       	= axes('Position', [clm row-.065 pl_dim(1) pl_dim(2)/5]); hold on
% ax0                       	= axes('Position', [clm row pl_dim]); hold on
% 
% % Extracte experimental data
% clear trg_acc trg_conf trg_coh trg_hit
% for iState = 1:size(t,1)
%     if t.trg_shown(iState) == false
%         continue
%     end
%     
%     for iTarget = 1:length(t.trg_ts{iState})
%         c                   = c+1;
%         trg_acc(c)          = t.trg_acc{iState}(iTarget);
%         trg_conf(c)         = t.trg_ecc{iState}(iTarget);
%         trg_coh(c)          = t.rdp_coh(iState);
%         trg_hit(c)          = t.trg_hit{iState}(iTarget);
%         
%         if trg_acc(c) < .5 && trg_hit(c) == 1
%             disp([num2str(iState) ' ' num2str(iTarget) ' ' num2str(c) ])
%         end
%     end
% end
% 
% % Calculate reward matrix
% acc                         = 0:.001:1;
% conf                        = 0:.001:1;
% rew                         = acc' .* conf;
% 
% % Determine arc width for each confidence level
% for j = 1:length(conf)
%     arc(j)                  = 180 - (180 * conf(j));
% end
% 
% % Cap arc width at target width (2dva == 12.7587deg at chosen position)
% aidx                        = arc < 12.7587;
% arc(aidx)                   = 12.7587;
% 
% % For each confidence level, calculate minimum accuracy required to hit
% % the target at given arc width - normalised values
% hit_width_acc               = 1 - ((arc/2) / 180);
% hit_width_acc(aidx)         = 1 - (12.7587/2)/180; % arc width fixed
% 
% % Remove position that cannot yield reward from reward matrix
% for iAcc = 1:length(acc)
%     indx                    = conf < hit_width_acc(iAcc);
%     rew(iAcc,indx)          = nan;
% end
% 
% % Plot reward matrix
% hold on
% im                          = imagesc(acc,conf,rew);
% ax0.XLabel.String           = 'Accuracy';
% ax0.YLabel.String           = 'Eccentricity';
% ax0.FontSize                = lb_fs;
% ax0.XLim                    = [0 1];
% ax0.YLim                    = [0 1];
% ax0.XLim                    = [0 1];
% ax0.XTick                   = [0:.2:1];
% ax0.YTick                   = [0:.2:1];
% ax0.XTickLabelRotation      = 0;
% ax0.YLabel.Position(1)      = -.4;
% ax0.XLabel.Position(2)      = -.22;
% 
% cb                          = colorbar;
% cb.Label.String             = 'Reward';
% cb.Location                 = 'eastoutside';
% colormap(cmap)
% cmap_coh                    = cool(size(snr,2));
% 
% for iCoh = 1:length(snr)
%     cidx                    = trg_coh == snr(iCoh);
%     sc(iCoh)                = scatter(trg_acc(cidx), trg_conf(cidx), 'filled');
%     sc(iCoh).CData        	= cmap_coh(iCoh,:);
%     sc(iCoh).SizeData      	= 15;
%     sc(iCoh).MarkerFaceAlpha = .9;
%     lg_str{iCoh}            = num2str(round(snr(iCoh)*100));
% end
% 
% ax0.Position                = [clm row pl_dim];
% 
% lg0                         = legend(sc,lg_str,'Location','northwest','NumColumns', 2);
% lg0.Box                     = 'off';
% lg0.Position(1:2)           = [.82 .61];
% 
% nBin                        = 40;
% axes(ax0h)
% [h, edg]                    = histcounts(trg_conf,nBin);
% cntr                        = edg(1:end-1) + diff(edg) ./ 2;
% st                          = stairs(-h,cntr);
% st.LineWidth                = lw/1.5;
% st.Color                    = [0 0 0];
% ax0h.YLim                   = [0 1];
% ax0h.XAxis.Visible          = 'off';
% ax0h.YAxis.Visible          = 'off';
% 
% axes(ax0v)
% [v, edg]                    = histcounts(trg_acc,nBin);
% cntr                        = edg(1:end-1) + diff(edg) ./ 2;
% st                          = stairs(cntr,-v);
% st.LineWidth                = lw/1.5;
% st.Color                    = [0 0 0];
% ax0v.XLim                   = [0 1];
% ax0v.XAxis.Visible          = 'off';
% ax0v.YAxis.Visible          = 'off';
% uistack(ax0v,'bottom')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accuracy-based psychometric curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% c = 0;
% for iSubj = 1:length(dyad_perf)
%     
%     if isempty(dyad_perf{iSubj})
%         continue
%     end
%     
%     if ~strcmp(dyad_perf{iSubj}.id, solo_perf{iSubj}.id) 
%         error('Subject confusion')
%     end
%     
%     c                     	= c+1;
%     crit                      = .875; % +/- 45 degrees
%     s                     	= cellfun(@(x) x > crit, solo_perf{iSubj}.acc_trg, 'UniformOutput', false);
%     d                     	= cellfun(@(x) x > crit, dyad_perf{iSubj}.acc_trg, 'UniformOutput', false);
%     
%     [param_solo(c,:) model_solo{c}] = fit_logistic_fct(solo_perf{iSubj}.carr, cellfun(@sum,s), cellfun(@length,s));
%     [param_dyad(c,:) model_dyad{c}] = fit_logistic_fct(dyad_perf{iSubj}.carr, cellfun(@sum,d), cellfun(@length,d));  
%     
% end
% 
% col                         = jet(size(param_solo,1));
% 
% ax18                       	= axes('Position', [colmn(2) height(5) .1 dim(2)]); hold on
% dat                         = param_solo(:,1) - param_dyad(:,1);
% xpos                        = randi([90 110],1,length(dat))./100;
% bx                          = boxplot(dat,'Color','k');
% sc                          = scatter(1:length(dat), dat, 'Marker','x');%, 'Color', col);
% ax18.YLim                   = [-.25 .25];
% ax18.XLim                   = [.8 1.2];
% ax18.XTick                  = [];
% ax18.YLabel.String          = 'd(Treshold)';
% ax18.FontSize               = lb_fs;
% 
% ax19                       	= axes('Position', [colmn(2)+.15 height(5) .1 dim(2)]); hold on
% bx                          = boxplot(param_solo(:,2) - param_dyad(:,2),'Color','k');
% ax19.XLim                   = [.8 1.2];
% ax19.XTick                  = [];
% ax19.YLabel.String          = 'd(Slope)';
% ax19.FontSize               = lb_fs;
% 
% ax20                       	= axes('Position', [colmn(1) height(5) dim]); hold on
% ex_id                       = 26;
% psolo                       = plot(model_solo{ex_id}(param_solo(ex_id,:),[0:.01:1]),'LineWidth',lw); hold on
% psolo.Color                 = [0 0 0 .5];
% pdyad                       = plot(model_dyad{ex_id}(param_dyad(ex_id,:),[0:.01:1]),'Color',[.75 0 0 .5],'LineWidth',lw);
% pdyad.Color                 = [.6 0 0 .5];
% ax20.XLabel.String          = 'Coherence';
% ax20.YLabel.String          = 'Hit rate [%]';
% ax20.FontSize               = lb_fs;
% lg                          = legend('solo','dyadic');
% lg.Location                 = 'southeast';
% lg.FontSize                 = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lag difference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dim                         = [.15 .15];
% ax6                         = axes('Position', [clmns(1) height(3) dim]); hold on
% 
% for iL = 1:length(lag)
%     pl(iL)               	= plot([1 2],[lag(iL,1) lag(iL,2)]);
%     pl(iL).Color          	= [.5 .5 .5 alp];
%     pl(iL).LineWidth      	= lw;
% end
% 
% bx                          = boxplot(lag, 'Colors', 'k');
% set(bx,'MarkerEdgeColor','k')
% set(bx, {'linew'},{lw})
% 
% ax6.YLabel.String           = 'Time [ms]';
% ax6.XTick                   = [1 2];
% ax6.XLim                    = [.5 2.5];
% ax6.XTickLabel              = {'Solo','Agnt'};
% ax6.FontSize                = lb_fs;
% ax6.XTickLabelRotation      = 0;
% ax6.Position                = [clmns(1) height(3) dim];
% ax6.Box                     = 'off';
% 
% function [paramsValues, PF] = fit_logistic_fct(snr, HitNo, OutOfNum)
% 
% PF = @PAL_Logistic;
% 
% %Threshold and Slope are free parameters, guess and lapse rate are fixed
% paramsFree = [1 1 0 0];  %1: free parameter, 0: fixed parameter
% 
% %Parameter grid defining parameter space through which to perform a
% %brute-force search for values to be used as initial guesses in iterative
% %parameter search.
% searchGrid.alpha = snr(1):.001:snr(end); % threshold
% searchGrid.beta = logspace(0,2,101); % slope
% searchGrid.gamma = 0;  % guess-rate
% searchGrid.lambda = 0.02;  % lapse-rate
% 
% %Perform fit
% disp('Fitting function.....');
% [paramsValues LL exitflag] = PAL_PFML_Fit(snr,HitNo,OutOfNum,searchGrid,paramsFree,PF);
% 
% disp('done:')
% message = sprintf('Threshold estimate: %6.4f',paramsValues(1));
% disp(message);
% message = sprintf('Slope estimate: %6.4f\r',paramsValues(2));
% disp(message);
% 
% %Number of simulations to perform to determine standard error
% B=400;
% disp('Determining standard errors.....');
% 
% [SD paramsSim LLSim converged] = PAL_PFML_BootstrapNonParametric(snr,HitNo, OutOfNum, [], paramsFree, B, PF,'searchGrid',searchGrid);
% 
% disp('done:');
% message = sprintf('Standard error of Threshold: %6.4f',SD(1));
% disp(message);
% message = sprintf('Standard error of Slope: %6.4f\r',SD(2));
% disp(message);
% 
% %Number of simulations to perform to determine Goodness-of-Fit
% B=1000;
% disp('Determining Goodness-of-fit.....');
% 
% [Dev pDev] = PAL_PFML_GoodnessOfFit(snr,HitNo,OutOfNum,paramsValues, paramsFree, B, PF, 'searchGrid', searchGrid);
% 
% disp('done:');
% message = sprintf('Deviance: %6.4f',Dev);
% disp(message);
% message = sprintf('p-value: %6.4f',pDev);
% disp(message);
% end

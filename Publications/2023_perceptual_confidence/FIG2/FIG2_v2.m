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

f                           = figure('units','normalized','position',[0 0 .5 1]);
height                      = [linspace(.82, .32,4) .1];
colmn                       = linspace(.075, .85,4);
lb_fs                       = 14;
lg_fs                       = 10;
lw                          = 3;
frme_ms                     = 1000/120;
alp                         = .25;
dim                         = [.18 .15];
col_dat                     = [0 0 0];
col_ci                      = [.3 0 0];
snr                         = solo_perf{end}.carr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: SOLO - Hit rate raw %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax1                         = axes('Position', [colmn(1) height(1) dim]); hold on
[ax1,pl]                    = plotData(ax1,solo_hir,true,lw,alp,col_dat,col_ci);
ax1.YLim                    = [10 70];
ax1.XLim                    = [1 size(solo_hir,2)];
ax1.XLabel.String           = 'Coherence [%]';
ax1.YLabel.String           = 'Hit rate [%]';
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
    for iCoh = 1:length(snr)
        clear cIdx
        cIdx             	= solo_cr{iSubj}.coh == snr(iCoh);
        solo_mlag(iSubj,iCoh) = mean(solo_cr{iSubj}.lag(cIdx));
    end
end

ax2                         = axes('Position', [colmn(1) height(4) dim]); hold on
[ax2,pl]                    = plotData(ax2,solo_mlag,false,lw,alp,col_dat,col_ci);
ax2.XLim                    = [1 size(solo_mlag,2)];
ax2.YLim                    = [300 1000];
ax2.XLabel.String           = 'Coherence [%]';
ax2.YLabel.String           = 'Lag [ms]';
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
ax3.YLabel.String           = 'Acc [%]';
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
ax4.YLabel.String           = 'Ecc [%]';
ax4.XTick                   = 1:length(snr);
ax4.XTickLabel              = round(snr,2)*100;
ax4.FontSize                = lb_fs;
ax4.XTickLabelRotation      = 0;
ax4.XAxis.Visible           = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: DYADIC - Hit rate raw %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax5                         = axes('Position', [colmn(2) height(1) dim]); hold on
[ax5,pl]                    = plotData(ax5,dyad_hir,true,lw,alp,col_dat,col_ci);
ax5.XLim                    = [1 size(dyad_hir,2)];
ax5.YLim                    = [10 70];
ax5.XLabel.String           = 'Coherence [%]';
ax5.YLabel.String           = [];
ax5.XTick                   = 1:length(snr);
ax5.XTickLabel              = round(snr,2)*100;
ax5.FontSize                = lb_fs;
ax5.XTickLabelRotation      = 0;
ax5.XAxis.Visible           = 'off';
ax5.YAxis.Visible           = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: Dyadic - Lag %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iSubj = 1:size(dyad_cr,2)
    if isempty(dyad_cr{iSubj})
        continue
    end
    
    for iCoh = 1:length(snr)
        clear cIdx
        cIdx             	= dyad_cr{iSubj}.coh == snr(iCoh);
        dyad_mlag(iSubj,iCoh) = mean(dyad_cr{iSubj}.lag(cIdx));
    end
end

ax6                         = axes('Position', [colmn(2) height(4) dim]); hold on
[ax6,pl]                    = plotData(ax6,dyad_mlag,false,lw,alp,col_dat,col_ci);
ax6.XLim                    = [1 size(dyad_mlag,2)];
ax6.YLim                    = [300 1000];
ax6.XLabel.String           = 'Coherence [%]';
ax6.YLabel.String           = 'Lag [ms]';
ax6.XTick                   = 1:length(snr);
ax6.XTickLabel              = round(snr,2)*100;
ax6.FontSize                = lb_fs;
ax6.XTickLabelRotation      = 0;
ax6.YAxis.Visible           = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: DYADIC - Avg accuracy raw %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax7                         = axes('Position', [colmn(2) height(2) dim]); hold on
[ax7,pl]                    = plotData(ax7,dyad_macc,true,lw,alp,col_dat,col_ci);
ax7.XLim                    = [1 size(dyad_macc,2)];
ax7.YLim                    = [30 100];
ax7.XLabel.String           = 'Coherence [%]';
ax7.YLabel.String           = [];
ax7.XTick                   = 1:length(snr);
ax7.XTickLabel              = round(snr,2)*100;
ax7.FontSize                = lb_fs;
ax7.XTickLabelRotation      = 0;
ax7.XAxis.Visible           = 'off';
ax7.YAxis.Visible           = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBPLOT: DYADIC - Avg eccentricity raw %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax8                         = axes('Position', [colmn(2) height(3) dim]); hold on
[ax8,pl]                    = plotData(ax8,dyad_mecc,true,lw,alp,col_dat,col_ci);
ax8.XLim                    = [1 size(dyad_mecc,2)];
ax8.YLim                    = [20 100];
ax8.XLabel.String           = 'Coherence [%]';
ax8.YLabel.String           = [];
ax8.XTick                   = 1:length(snr);
ax8.XTickLabel              = round(snr,2)*100;
ax8.FontSize                = lb_fs;
ax8.XTickLabelRotation      = 0;
ax8.XAxis.Visible           = 'off';
ax8.YAxis.Visible           = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBPLOT: Hit rate comparison %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cnt = 0;
for iSubj = 1:30%size(sbj_lst,1)
    if ~isempty(dyad_perf{iSubj})
        cnt = cnt+1;
        df_hir(cnt,:) = solo_perf{iSubj}.hir - dyad_perf{iSubj}.hir;
        df_hir_pool(cnt,:) = solo_perf{iSubj}.hir_pool - dyad_perf{iSubj}.hir_pool;
    end
end

ax13                      	= axes('Position', [colmn(3) height(1) dim]); hold on
pt                        	= patch([0 8 8 0], [0 0 25 25], [.75 1 .75], 'FaceAlpha',.1, 'EdgeColor','none');
pt                       	= patch([0 8 8 0], [-25 -25 0 0], [1 .75 .75], 'FaceAlpha',.1, 'EdgeColor','none');
tx                        	= text(1.25,18, 'Dyad > Solo', 'FontSize', lb_fs, 'Color', 'k');
tx                       	= text(1.25,-18, 'Dyad < Solo', 'FontSize', lb_fs, 'Color', 'k');

[ax13,pl]                   = plotData(ax13,df_hir,true,lw,alp,col_dat,col_ci);
lm                          = line([1 7],[0 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth',lw);
ax13.XLim                   = [1 size(df_hir,2)];
ax13.FontSize               = lb_fs;
ax13.YLabel.String          = 'Difference';
ax13.XAxis.Visible          = 'off';

% ax17                       	= axes('Position', [colmn(4) height(1) .1 dim(2)]); hold on
% ax17                      	= plotBar(ax17, df_hir, lb_fs,snr);

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

for iSubj = 1:size(dyad_perf,2)
    
    d.sp                  	= solo_perf{iSubj};                             % Performance data
    d.sc                  	= solo_cr{iSubj};                               % Correlation data
    
    sIdx                    = cellfun(@(x) strcmp(x,solo_perf{iSubj}.id),ply_id);
    
    if sum(sIdx) == 0
        continue
    end
    
    d.dp                  	= dyad_perf{sIdx};
    d.dc                  	= dyad_cr{iSubj};
         
    if isempty(d.dc)
        continue
    end
    
    for iCoh = 1:length(snr)
        auc_acc(iSubj,iCoh)  	= f_auroc(d.sp.acc_trg{iCoh},d.dp.acc_trg{iCoh});
%         auc_acc(iSubj,iCoh)  	= f_auroc(d.sp.acc_state{iCoh},d.dp.acc_state{iCoh});
        auc_ecc(iSubj,iCoh)   	= f_auroc(d.sp.ecc{iCoh},d.dp.ecc{iCoh});
        auc_score(iSubj,iCoh)  	= f_auroc(d.sp.trg_score{iCoh},d.dp.trg_score{iCoh});
        
        auc_cc(iSubj,iCoh)    	= f_auroc(d.sc.cc(snr(iCoh) == d.sc.coh),d.dc.cc(snr(iCoh) == d.dc.coh));
        auc_xcp(iSubj,iCoh)   	= f_auroc(d.sc.posPk(snr(iCoh) == d.sc.coh),d.dc.posPk(snr(iCoh) == d.dc.coh));
        auc_xc(iSubj,iCoh)    	= f_auroc(d.sc.maxR(snr(iCoh) == d.sc.coh),d.dc.maxR(snr(iCoh) == d.dc.coh));
        auc_lag(iSubj,iCoh)    	= f_auroc(d.sc.lag(snr(iCoh) == d.sc.coh),d.dc.lag(snr(iCoh) == d.dc.coh));
        
        kk(iSubj,iCoh) = mean(d.sc.lag(snr(iCoh) == d.sc.coh))-mean(d.dc.lag(snr(iCoh) == d.dc.coh));
        
        p_acc(iSubj,iCoh)      	= ranksum(d.sp.acc_trg{iCoh},d.dp.acc_trg{iCoh});
%         p_acc(iSubj,iCoh)      	= ranksum(d.sp.acc_state{iCoh},d.dp.acc_state{iCoh});
        p_ecc(iSubj,iCoh)      	= ranksum(d.sp.ecc{iCoh},d.dp.ecc{iCoh});
        p_score(iSubj,iCoh)  	= ranksum(d.sp.trg_score{iCoh},d.dp.trg_score{iCoh});
        p_cc(iSubj,iCoh)      	= ranksum(d.sc.cc(snr(iCoh) == d.sc.coh),d.dc.cc(snr(iCoh) == d.dc.coh));
        p_xcp(iSubj,iCoh)      	= ranksum(d.sc.posPk(snr(iCoh) == d.sc.coh),d.dc.posPk(snr(iCoh) == d.dc.coh));
        p_xc(iSubj,iCoh)       	= ranksum(d.sc.maxR(snr(iCoh) == d.sc.coh),d.dc.maxR(snr(iCoh) == d.dc.coh));
        p_lag(iSubj,iCoh)      	= ranksum(d.sc.lag(snr(iCoh) == d.sc.coh),d.dc.lag(snr(iCoh) == d.dc.coh));
    end
end

ax9                       	= axes('Position', [colmn(3) height(4) dim]); hold on
pt                         	= patch([0 8 8 0], [.5 .5 1 1], [.75 1 .75], 'FaceAlpha',.1, 'EdgeColor','none');
pt                         	= patch([0 8 8 0], [0 0 .5 .5], [1 .75 .75], 'FaceAlpha',.1, 'EdgeColor','none');

ax10                     	= axes('Position', [colmn(3) height(2) dim]); hold on
pt                         	= patch([0 8 8 0], [.5 .5 1 1], [.75 1 .75], 'FaceAlpha',.1, 'EdgeColor','none');
pt                         	= patch([0 8 8 0], [0 0 .5 .5], [1 .75 .75], 'FaceAlpha',.1, 'EdgeColor','none');
tx                        	= text(1.25,.9, 'Dyad > Solo', 'FontSize', lb_fs, 'Color', 'k');
tx                       	= text(1.25,.1, 'Dyad < Solo', 'FontSize', lb_fs, 'Color', 'k');

ax11                       	= axes('Position', [colmn(3) height(3) dim]); hold on
pt                         	= patch([0 8 8 0], [.5 .5 1 1], [.75 1 .75], 'FaceAlpha',.1, 'EdgeColor','none');
pt                         	= patch([0 8 8 0], [0 0 .5 .5], [1 .75 .75], 'FaceAlpha',.1, 'EdgeColor','none');

ax9                     	= plotAUROC(ax9,auc_lag,'AUROC',lb_fs,snr,alp,lw,col_dat,col_ci);
ax9.XAxis.Visible         	= 'on';
ax10                     	= plotAUROC(ax10,auc_acc,'AUROC',lb_fs,snr,alp,lw,col_dat,col_ci);
ax11                    	= plotAUROC(ax11,auc_ecc,'AUROC',lb_fs,snr,alp,lw,col_dat,col_ci);


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Accuracy-based psychometric curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 0;
for iSubj = 1:length(dyad_perf)
    
    if isempty(dyad_perf{iSubj})
        continue
    end
    
    if ~strcmp(dyad_perf{iSubj}.id, solo_perf{iSubj}.id) 
        error('Subject confusion')
    end
    
    c                     	= c+1;
    crit                      = .875; % +/- 45 degrees
    s                     	= cellfun(@(x) x > crit, solo_perf{iSubj}.acc_trg, 'UniformOutput', false);
    d                     	= cellfun(@(x) x > crit, dyad_perf{iSubj}.acc_trg, 'UniformOutput', false);
    
    [param_solo(c,:) model_solo{c}] = fit_logistic_fct(solo_perf{iSubj}.carr, cellfun(@sum,s), cellfun(@length,s));
    [param_dyad(c,:) model_dyad{c}] = fit_logistic_fct(dyad_perf{iSubj}.carr, cellfun(@sum,d), cellfun(@length,d));  
    
end

col                         = jet(size(param_solo,1));

ax18                       	= axes('Position', [colmn(2) height(5) .1 dim(2)]); hold on
dat                         = param_solo(:,1) - param_dyad(:,1);
xpos                        = randi([90 110],1,length(dat))./100;
bx                          = boxplot(dat,'Color','k');
sc                          = scatter(1:length(dat), dat, Color, col, Marker,'x');
ax18.YLim                   = [-.25 .25];
ax18.XLim                   = [.8 1.2];
ax18.XTick                  = [];
ax18.YLabel.String          = 'd(Treshold)';
ax18.FontSize               = lb_fs;

ax19                       	= axes('Position', [colmn(2)+.15 height(5) .1 dim(2)]); hold on
bx                          = boxplot(param_solo(:,2) - param_dyad(:,2),'Color','k');
ax19.XLim                   = [.8 1.2];
ax19.XTick                  = [];
ax19.YLabel.String          = 'd(Slope)';
ax19.FontSize               = lb_fs;

ax20                       	= axes('Position', [colmn(1) height(5) dim]); hold on
ex_id                       = 26;
psolo                       = plot(model_solo{ex_id}(param_solo(ex_id,:),[0:.01:1]),'LineWidth',lw); hold on
psolo.Color                 = [0 0 0 .5];
pdyad                       = plot(model_dyad{ex_id}(param_dyad(ex_id,:),[0:.01:1]),'Color',[.75 0 0 .5],'LineWidth',lw);
pdyad.Color                 = [.6 0 0 .5];
ax20.XLabel.String          = 'Coherence';
ax20.YLabel.String          = 'Hit rate [%]';
ax20.FontSize               = lb_fs;
lg                          = legend('solo','dyadic');
lg.Location                 = 'southeast';
lg.FontSize                 = 10;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Annotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dest_dir                    = '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG2/';
ax0                         = axes('Position',[0 0 1 1],'Visible','off');

text(0.12,.98, 'Solo', 'Parent', ax0, 'FontSize', 22, 'Color', 'k')
text(0.37,.98, 'Dyadic', 'Parent', ax0, 'FontSize', 22, 'Color', 'k')
text(0.61,.98, 'Contrast', 'Parent', ax0, 'FontSize', 22, 'Color', 'k')
text(0.85,.98, 'Stats', 'Parent', ax0, 'FontSize', 22, 'Color', 'k')

print(f, [dest_dir '/FIG2'], '-r500', '-dpng');

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
clear tscore score_cum score_hi score_coh trg_states
for iState = 1:size(in.trg_ts,1)
    for iTrg = 1:length(in.trg_ts{iState})
        try
            c                   = c +1;
            score_hi(c)         = in.trg_hit{iState}(iTrg);
            score(c)            = double(in.trg_score{iState}(iTrg));
            score_coh(c)        = in.rdp_coh(iState);
        catch
            warning(['Skipped state/target: ' num2str(iState) '/' num2str(iTrg)])
        end
    end
end

trg_states                 	= in.trg_hit(logical(in.trg_shown));
out.hir_pool                = sum(cellfun(@sum,trg_states)) / sum(cellfun(@numel,trg_states));

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
    pl(iL)                  = plot(dat(iL,:), 'LineWidth', lw/3, 'Color', [.5 .5 .5 alp]);
end

% Overlay confidence intervals
fl                          = fill(x_spacing,ci,col_ci,'EdgeColor','none', 'FaceAlpha', alp);

% Plot mean curve
pl                          = plot(mean(dat), 'LineWidth', lw, 'Color', col_dat);
end

function ax = plotAUROC(ax,dat,str,lb_fs,snr,alp,lw,col_dat,col_ci)

axes(ax); hold on

% Remoe zero rows
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
ax.YLim                  	= [0 1];
ax.FontSize               	= lb_fs;
ax.YLabel.String           	= str;
ax.XLabel.String           	= 'Coherence [%]';
ax.XTick                 	= 1:7;
ax.YTick                  	= [.2 .5 .8];
ax.XTickLabel            	= round(snr,2).*100;
ax.XTickLabelRotation      	= 0;
ax.XAxis.Visible            = 'off';

n                           = .05;
crit                        = [.5-n .5+n];
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
bp.FaceColor                = [0 .6 0];

bn                          = bar( -(sum(neg)/all) .*100);
bn.FaceColor                = [.6 0 0];

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

function [paramsValues, PF] = fit_logistic_fct(snr, HitNo, OutOfNum)

PF = @PAL_Logistic;

%Threshold and Slope are free parameters, guess and lapse rate are fixed
paramsFree = [1 1 0 0];  %1: free parameter, 0: fixed parameter

%Parameter grid defining parameter space through which to perform a
%brute-force search for values to be used as initial guesses in iterative
%parameter search.
searchGrid.alpha = snr(1):.001:snr(end); % threshold
searchGrid.beta = logspace(0,2,101); % slope
searchGrid.gamma = 0;  % guess-rate
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
end
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

nSample                     = 29;                                           % Time window size [samples]
nLag                        = 150;                                          % Cross-correlation lag
s_cnt                       = 0;                                            % Subject counter
d_cnt                       = 0;                                            % Dyad counter

%% Extract data

% For all subjects
for iSubj = 1:length(sbj_lst)
    
    disp(['Processing subject: ' sbj_lst{iSubj}])
    
    data_pth                     	= [pth sbj_lst{iSubj} '/summary/'];   	% Data path
    solo_tbl                     	= [];                                   % Initialise
    agnt_tbl                     	= [];                                   
    dyad_tbl                     	= [];
    dyad_pc_tbl                     = [];
    tmp_solo                        = [];
    tmp_agnt                        = [];
    tmp_dyad                        = [];
    tmp_dyad_pc                     = [];
    sc_cnt                       	= 0;                                    % Reset cycle counter
    ac_cnt                       	= 0;                                   
    dc_cnt                       	= 0;
    dpc_cnt                       	= 0;
    
    if isdir(data_pth)
        cd(data_pth)
        mat_files               	= dir('*.mat');                         % Get all .mat files in directory
        
        % For all files in directory
        for iFile = 1:length(mat_files)

            % Solo experiments
            if contains(mat_files(iFile).name,'CPRsolo')
                % Load pre-processed data table
                tmp_solo             = load(mat_files(iFile).name);
                
                % Exclude sessions with different coherence level
                if sum(unique(tmp_solo.t.rdp_coh) < .1) > 2
                    continue
                end
                
                % Concatenate session summary tables
                solo_tbl           	= [solo_tbl; tmp_solo.t];                % Organise solo data
            end
            
            % Computer - human dyads
            if contains(mat_files(iFile).name,'CPRagent') && ~contains(mat_files(iFile).name,'agnt') 
                % Load pre-processed data table
                tmp_dyad_pc      	= load(mat_files(iFile).name);
                
                % Concatenate session summary tables
                dyad_pc_tbl       	= [dyad_pc_tbl; tmp_dyad_pc.t];       	% Organise dyad[human-computer] data
            end
            
            % Computer player performance
            if contains(mat_files(iFile).name,'agnt_CPRagent')
                % Load pre-processed data table
                tmp_agnt             = load(mat_files(iFile).name);

                % Concatenate session summary tables
                agnt_tbl           	= [agnt_tbl; tmp_agnt.t];                % Organise computer player data
            end
        end
        
        % For human - human dyads
        dyad_cnt                    = 0;
        for iDyad = 19:71
            
            if iDyad == 36 || iDyad == 55
                continue
            end
            
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
                    dyad_tbl      	= [dyad_tbl; tmp_tbl.t];                % Organise dyad[human-human] data
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% PERFORMANCE AND CORRELATION ANALYSIS %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if ~isempty(solo_tbl)
            %% Performance [time window analysis]
            solo_perf{iSubj}      	= response_readout(solo_tbl, nSample);

            % Extract stimulus cycles for correlation analysis
            cyc_boundary        	= [0; find(diff(solo_tbl.cyc_no) ~= 0)];
            [tmp_solo, sc_cnt]   	= extractCycles(cyc_boundary, solo_tbl, tmp_solo, sc_cnt);
            
            % Correlation analysis
            [solo_cr{iSubj},ps]     = CPR_correlation_analysis_WIP(tmp_solo, nLag, false);
            
            solo_perf{iSubj}.id     = lower(sbj_lst{iSubj});
            solo_cr{iSubj}.id       = lower(sbj_lst{iSubj});
        end
        
        if ~isempty(agnt_tbl)
            % Performance [time window analysis]
            agnt_perf{iSubj}      	= response_readout(agnt_tbl, nSample);

            % Extract stimulus cycles for correlation analysis
            cyc_boundary        	= [0; find(diff(agnt_tbl.cyc_no) ~= 0)];
            [tmp_agnt, ac_cnt]   	= extractCycles(cyc_boundary, agnt_tbl, tmp_agnt, ac_cnt);
            
            % Correlation analysis
            [agnt_cr{iSubj},ps]     = CPR_correlation_analysis_WIP(tmp_agnt, nLag, false);
            
            agnt_perf{iSubj}.id     = lower(sbj_lst{iSubj});
            agnt_cr{iSubj}.id       = lower(sbj_lst{iSubj});
        end
        
        if ~isempty(dyad_pc_tbl)
            % Performance [time window analysis]
            dyad_pc_perf{iSubj}   	= response_readout(dyad_pc_tbl, nSample);
            
            % Extract stimulus cycles for correlation analysis
            cyc_boundary          	= [0; find(diff(dyad_pc_tbl.cyc_no) ~= 0)];
            [tmp_dyad_pc, dpc_cnt]	= extractCycles(cyc_boundary, dyad_pc_tbl, tmp_dyad_pc, dpc_cnt);
            
            % Correlation analysis
            [dyad_pc_cr{iSubj},ps]	= CPR_correlation_analysis_WIP(tmp_dyad_pc, nLag, false);
            
            dyad_pc_perf{iSubj}.id 	= lower(sbj_lst{iSubj});
            dyad_pc_cr{iSubj}.id  	= lower(sbj_lst{iSubj});
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

%%% Save variables for plotting %%%
dest_dir = ['/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/'];
mkdir(dest_dir);

save([dest_dir 'solo_performance.mat'], 'solo_perf', '-v7.3');
save([dest_dir 'comp_performance.mat'], 'agnt_perf', '-v7.3');
save([dest_dir 'dyad_human_comp_performance.mat'], 'dyad_pc_perf', '-v7.3');
save([dest_dir 'dyad_human_human_performance.mat'], 'dyad_perf', '-v7.3');
save([dest_dir 'solo_correlation.mat'], 'solo_cr', '-v7.3');
save([dest_dir 'comp_correlation.mat'], 'agnt_cr', '-v7.3');
save([dest_dir 'dyad_human_human_correlation.mat'], 'dyad_cr', '-v7.3');
save([dest_dir 'dyad_human_comp_correlation.mat'], 'dyad_pc_cr', '-v7.3');

%% Import dyadic data [human-human]

disp('DYAD_WISE ANALYSIS')

pth                         = '/Volumes/T7_Shield/CPR_psychophysics/';      % Local hard drive
dyad_cnt                    = 0;
nSample                     = 30;
nLag                        = 150;

for iDyad = 19:71
    
    disp(['Dyad:' num2str(iDyad)])
    
    if iDyad == 36 || iDyad == 55
        continue
    end
            
    cd([pth ['Dyad' num2str(iDyad)] '/summary/'])
    mat_files              	= dir('*.mat');
    
    if isempty(mat_files)
        continue
    end
    
    dyad_cnt                = dyad_cnt+1;
    tmp1                    = load(mat_files(1).name);
    tmp3                    = load(mat_files(3).name);
    id_dyad(dyad_cnt,:)    	= [{tmp1.t.ID(1)}, {tmp3.t.ID(1)}];
    n_dyad(dyad_cnt,:)      = iDyad;
    cnt                     = 0;
    
    for iPly = 1:2
        if iPly == 1
            tbl             = tmp1.t;
        else
            tbl             = tmp3.t;
        end
                       
        % Performance [time window analysis]
        dyad_pw_perf{dyad_cnt,iPly}     = response_readout(tbl, nSample);
        
        % Extract stimulus cycles for correlation analysis
        cyc_boundary                    = [0; find(diff(tbl.cyc_no) ~= 0)];
        [tmp_dyad, cnt]                 = extractCycles(cyc_boundary, tbl, [], cnt);

        % Correlation analysis
        [dyad_pw_cr{dyad_cnt,iPly},ps] 	= CPR_correlation_analysis_WIP(tmp_dyad, nLag, false);
        
        dyad_pw_perf{dyad_cnt,iPly}.id  = lower(tbl.ID(1));
        dyad_pw_cr{dyad_cnt,iPly}.id    = lower(tbl.ID(1));   
    end 
end

save([dest_dir 'dyad_pairwise_performance.mat'], 'dyad_pw_perf', '-v7.3');
save([dest_dir 'dyad_pairwise_correlation.mat'], 'dyad_pw_cr', '-v7.3');

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

dte                         = unique(score_dte);
dte(ismissing(dte))         = [];

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
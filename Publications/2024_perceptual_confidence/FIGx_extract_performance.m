% Add relevant directories
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/preprocessing/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/mat_to_summary/
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab
addpath /Users/fschneider/Documents/MATLAB/cbrewer/

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
    comp_tbl                     	= [];                                   
    dyad_tbl                     	= [];
    solo_tbl_filt                   = [];                                   % Initialise
    comp_tbl_filt                	= [];
    dyad_tbl_filt                	= [];
    hc_dyad_tbl_filt              	= [];
    hc_dyad_tbl                     = [];
    tmp_solo                        = [];
    tmp_comp                        = [];
    tmp_dyad                        = [];
    tmp_hc_dyad                     = [];
    sc_cnt                       	= 0;                                    % Reset cycle counter
    ac_cnt                       	= 0;                                   
    dc_cnt                       	= 0;
    dhc_cnt                       	= 0;
    
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
                solo_tbl_filt       = [solo_tbl; tmp_solo.t(~logical(tmp_solo.t.fix_flag),:)]; % filtered for fixation breaks
                solo_tbl           	= [solo_tbl; tmp_solo.t];
            end
            
            % Computer - human dyads
            if contains(mat_files(iFile).name,'CPRagent') && ~contains(mat_files(iFile).name,'agnt') 
                % Load pre-processed data table
                tmp_hc_dyad      	= load(mat_files(iFile).name);
                
                % Concatenate session summary tables
                hc_dyad_tbl_filt   	= [hc_dyad_tbl_filt; tmp_hc_dyad.t(~logical(tmp_hc_dyad.t.fix_flag),:)]; % Organise dyad [human-computer] data
                hc_dyad_tbl       	= [hc_dyad_tbl; tmp_hc_dyad.t];
            end
            
            % Computer player performance
            if contains(mat_files(iFile).name,'agnt_CPRagent')
                % Load pre-processed data table
                tmp_comp             = load(mat_files(iFile).name);

                % Concatenate session summary tables
                comp_tbl           	 = [comp_tbl; tmp_comp.t];                % Organise computer player data
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
            tmp3                    = split(mat_files(end).name,'_');
            id_dyad(dyad_cnt,:)    	= [{tmp1{2}}, {tmp3{2}}];
            n_dyad(dyad_cnt,:)      = iDyad;
            
            % Check if subject part of this dyad
            for iFile = 1:size(mat_files,1)
                is_player(iFile) = contains(mat_files(iFile).name,lower(sbj_lst{iSubj}));
                
                if is_player(iFile)
                    tmp_tbl      	= load(mat_files(iFile).name);
                    dyad_tbl_filt   = [dyad_tbl_filt; tmp_tbl.t(~logical(tmp_tbl.t.fix_flag),:)]; % Organise dyad[human-human] data
                    dyad_tbl      	= [dyad_tbl; tmp_tbl.t];               
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% PERFORMANCE AND CORRELATION ANALYSIS %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                             
%         % Visualise target eccentricity
%         plot_trg_response(solo_tbl)
%         plot_trg_response(dyad_tbl)
%         plot_trg_response(hc_dyad_tbl)
        
        if ~isempty(solo_tbl)
            %%% Performance [time window analysis]
            solo_perf{iSubj}      	= response_readout(solo_tbl_filt, nSample);

            % Extract stimulus cycles for correlation analysis
            cyc_boundary        	= [0; find(diff(solo_tbl.cyc_no) ~= 0)];
            [tmp_solo, sc_cnt]   	= extractCycles(cyc_boundary, solo_tbl, tmp_solo, sc_cnt);
            
            %%% Correlation analysis
            [solo_cr{iSubj},ps]     = CPR_correlation_analysis(tmp_solo, nLag, false);
            
            % Add ID
            solo_perf{iSubj}.id     = lower(sbj_lst{iSubj});
            solo_cr{iSubj}.id       = lower(sbj_lst{iSubj});
        end
        
        if ~isempty(comp_tbl)
            %%% Performance [time window analysis]
            comp_perf{iSubj}      	= response_readout(comp_tbl, nSample);

            % Extract stimulus cycles for correlation analysis
            cyc_boundary        	= [0; find(diff(comp_tbl.cyc_no) ~= 0)];
            [tmp_comp, ac_cnt]   	= extractCycles(cyc_boundary, comp_tbl, tmp_comp, ac_cnt);
            
            %%% Correlation analysis
            [comp_cr{iSubj},ps]     = CPR_correlation_analysis(tmp_comp, nLag, false);
            
            % Add ID
            comp_perf{iSubj}.id     = lower(sbj_lst{iSubj});
            comp_cr{iSubj}.id       = lower(sbj_lst{iSubj});
        end
        
        if ~isempty(hc_dyad_tbl)
            %%% Performance [time window analysis]
            hc_dyad_perf{iSubj}   	= response_readout(hc_dyad_tbl_filt, nSample);
            
            % Extract stimulus cycles for correlation analysis
            cyc_boundary          	= [0; find(diff(hc_dyad_tbl.cyc_no) ~= 0)];
            [tmp_hc_dyad, dhc_cnt]	= extractCycles(cyc_boundary, hc_dyad_tbl, tmp_hc_dyad, dhc_cnt);
            
            %%% Correlation analysis
            [hc_dyad_cr{iSubj},ps]	= CPR_correlation_analysis(tmp_hc_dyad, nLag, false);
                       
            % Add ID
            hc_dyad_perf{iSubj}.id 	= lower(sbj_lst{iSubj});
            hc_dyad_cr{iSubj}.id  	= lower(sbj_lst{iSubj}); 
        end
        
        if ~isempty(dyad_tbl)
            %%% Performance [time window analysis]
            dyad_perf{iSubj}      	= response_readout(dyad_tbl_filt, nSample);
            
            % Extract stimulus cycles for correlation analysis
            cyc_boundary        	= [0; find(diff(dyad_tbl.cyc_no) ~= 0)];
            [tmp_dyad, dc_cnt]   	= extractCycles(cyc_boundary, dyad_tbl, tmp_dyad, dc_cnt);
            
            %%% Correlation analysis
            [dyad_cr{iSubj},ps]     = CPR_correlation_analysis(tmp_dyad, nLag, false);
            
            % Add ID
            dyad_perf{iSubj}.id     = lower(sbj_lst{iSubj});
            dyad_cr{iSubj}.id       = lower(sbj_lst{iSubj});
        end
    end
end

%%% Save variables for plotting %%%
% dest_dir = ['/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/'];
dest_dir = ['/Users/fschneider/ownCloud/var_plot/'];
mkdir(dest_dir);


save([dest_dir 'solo_performance.mat'], 'solo_perf', '-v7.3');
save([dest_dir 'comp_performance.mat'], 'comp_perf', '-v7.3');
save([dest_dir 'hc_dyad_performance.mat'], 'hc_dyad_perf', '-v7.3');
save([dest_dir 'hh_dyad_performance.mat'], 'dyad_perf', '-v7.3');
save([dest_dir 'solo_correlation.mat'], 'solo_cr', '-v7.3');
save([dest_dir 'comp_correlation.mat'], 'comp_cr', '-v7.3');
save([dest_dir 'hh_dyad_correlation.mat'], 'dyad_cr', '-v7.3');
save([dest_dir 'hc_dyad_correlation.mat'], 'hc_dyad_cr', '-v7.3');

%% Import dyadic data [human-human]

disp('HH-DYAD_WISE ANALYSIS')

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
    
    % Extract player identity
    fname_split             = split(mat_files(1).name,'_');
    p1_str                 	= fname_split{2};
    fname_split             = split(mat_files(end).name,'_');
    p2_str                 	= fname_split{2};
    dyad_cnt                = dyad_cnt+1;
    id_dyad(dyad_cnt,:)    	= [{p1_str}, {p2_str}];
    n_dyad(dyad_cnt,:)      = iDyad;
    
    % Initialise variables
    fcnt                    = 0;
    cnt                     = 0;
    p1_flag                 = [];
    dyad_file               = [];
    hum1_tbl                = [];
    hum1_tbl_filt           = [];
    hum2_tbl                = [];
    hum2_tbl_filt           = [];
        
    % Get all dyadic files
    for iFile = 1:length(mat_files)
        if contains(mat_files(iFile).name,'CPRdyadic')
            fcnt                = fcnt+1;
            dyad_file{fcnt}     = mat_files(iFile).name;
            if contains(mat_files(iFile).name,p1_str)
                p1_flag(fcnt)   = true;
            else
                p1_flag(fcnt)   = false;
            end
        end
    end
    
    if isempty(dyad_file)
        continue
    end

    % Concatenate data tables for each player
    for iHuman = find(p1_flag)
        tmp                 = load(dyad_file{iHuman});
        hum1_tbl_filt        = [hum1_tbl_filt; tmp.t(~logical(tmp.t.fix_flag),:)];
        hum1_tbl             = [hum1_tbl; tmp.t];
    end
    
    for iHuman = find(~p1_flag)
        tmp                 = load(dyad_file{iHuman});
        hum2_tbl_filt        = [hum2_tbl_filt; tmp.t(~logical(tmp.t.fix_flag),:)];
        hum2_tbl             = [hum2_tbl; tmp.t];
    end
    
    % Analyse behavior of each player
    for iPly = 1:2
        if iPly == 1
            tbl             = hum1_tbl;
            tbl_filt        = hum1_tbl_filt;
        else
            tbl             = hum2_tbl;
            tbl_filt        = hum2_tbl_filt;
        end
                       
        % Performance [time window analysis]
        dyad_pw_perf{dyad_cnt,iPly}     = response_readout(tbl_filt, nSample);
        
        % Extract stimulus cycles for correlation analysis
        cyc_boundary                    = [0; find(diff(tbl.cyc_no) ~= 0)];
        [tmp_dyad, cnt]                 = extractCycles(cyc_boundary, tbl, [], cnt);

        % Correlation analysis
        [dyad_pw_cr{dyad_cnt,iPly},ps] 	= CPR_correlation_analysis(tmp_dyad, nLag, false);
        
        dyad_pw_perf{dyad_cnt,iPly}.id  = lower(tbl.ID(1));
        dyad_pw_cr{dyad_cnt,iPly}.id    = lower(tbl.ID(1));   
    end 
end

save([dest_dir 'hh_dyad_pairwise_performance.mat'], 'dyad_pw_perf', '-v7.3');
save([dest_dir 'hh_dyad_pairwise_correlation.mat'], 'dyad_pw_cr', '-v7.3');

%% Import dyadic data [computer-human]

disp('HC-DYAD_WISE ANALYSIS')

pth                         = '/Volumes/T7_Shield/CPR_psychophysics/';      % Local hard drive
dyad_cnt                    = 0;
nSample                     = 30;
nLag                        = 150;

for iSubj = 1:length(sbj_lst)
    
    disp(['iSubj:' num2str(iSubj)])
   
    cd([pth [sbj_lst{iSubj}] '/summary/'])
    mat_files              	= dir('*.mat');
    
    cnt                     = 0;
    hc_dyad_file            = [];
    comp_tbl                = [];
    hum_tbl                 = [];
    hum_tbl_filt         	= [];
    comp_flag               = [];
    
    for iFile = 1:length(mat_files)
        if contains(mat_files(iFile).name,'CPRagent') 
            cnt                 = cnt+1;
            hc_dyad_file{cnt}   = mat_files(iFile).name;
            if contains(mat_files(iFile).name,'agnt')
                comp_flag(cnt) = true;
            else
                comp_flag(cnt) = false;
            end
        end
    end
    
    if isempty(hc_dyad_file)
        continue
    end
        
    for iComputer = find(comp_flag)
        tmp                 = load(hc_dyad_file{iComputer});
        comp_tbl            = [comp_tbl; tmp.t];
    end
    
    for iHuman = find(~comp_flag)
        tmp                 = load(hc_dyad_file{iHuman});
        hum_tbl_filt        = [hum_tbl_filt; tmp.t(~logical(tmp.t.fix_flag),:)];
        hum_tbl             = [hum_tbl; tmp.t];
    end
    
    dyad_cnt                = dyad_cnt+1;
    id_dyad(dyad_cnt,:)    	= [{comp_tbl.ID(end)}, {hum_tbl.ID(end)}];
    cnt                     = 0;

    for iPly = 1:2
        if iPly == 1
            tbl             = comp_tbl;
            tbl_filt        = comp_tbl;
        else
            tbl             = hum_tbl;
            tbl_filt        = hum_tbl_filt;
        end
                       
        % Performance [time window analysis]
        hc_dyad_pw_perf{dyad_cnt,iPly}     = response_readout(tbl_filt, nSample);
        
        % Extract stimulus cycles for correlation analysis
        cyc_boundary                        = [0; find(diff(tbl.cyc_no) ~= 0)];
        [tmp_dyad, cnt]                     = extractCycles(cyc_boundary, tbl, [], cnt);

        % Correlation analysis
        [hc_dyad_pw_cr{dyad_cnt,iPly},ps] 	= CPR_correlation_analysis(tmp_dyad, nLag, false);
        
        hc_dyad_pw_perf{dyad_cnt,iPly}.id   = lower(tbl.ID(1));
        hc_dyad_pw_cr{dyad_cnt,iPly}.id     = lower(tbl.ID(1));   
    end 
end

save([dest_dir 'hc_dyad_pairwise_performance.mat'], 'hc_dyad_pw_perf', '-v7.3');
save([dest_dir 'hc_dyad_pairwise_correlation.mat'], 'hc_dyad_pw_cr', '-v7.3');

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


%%% For all dates %%%
dte                         = unique(score_dte);
dte(ismissing(dte))         = [];

for iDate = 1:length(dte)
    clear dte_idx score_cum dte_score
    dte_idx                 = score_dte == dte(iDate);
    out.score_final_exp(iDate) = unique(score_exp(dte_idx));
    dte_score               = score(dte_idx);
    out.score_final(iDate)  = sum(dte_score(~isnan(dte_score)));
    out.score_norm(iDate)   = out.score_final(iDate) ./ length(dte_score(~isnan(dte_score)));
end

%%% Target-wise %%%
in.trg_score             	= cellfun(@double,in.trg_score,'UniformOutput', false);
out.trg_all_outc          	= logical(cell2mat(in.trg_hit(logical(in.trg_shown))'));
out.trg_all_score          	= cell2mat(in.trg_score(logical(in.trg_shown))');
out.trg_all_ecc          	= cell2mat(in.trg_ecc(logical(in.trg_shown))');
out.trg_all_acc           	= cell2mat(in.trg_acc(logical(in.trg_shown))');

trg_n                       = cellfun(@length,in.trg_ecc(logical(in.trg_shown)));
trg_coh                     = in.rdp_coh(logical(in.trg_shown));
tmp_trg_ts                  = in.trg_ts(logical(in.trg_shown));
tmp_frame1_ts             	= cellfun(@(x) x(1), in.frme_ts(logical(in.trg_shown)));
    
for j = 1:length(trg_n)
    tmp_trg_all_coh{j}   	= repmat(trg_coh(j),1,trg_n(j));
    tmp_trg_ts_state{j} 	= double((tmp_trg_ts{j} - tmp_frame1_ts(j)) ./ 1e3);  
end
out.trg_all_coh            	= cell2mat(tmp_trg_all_coh);
out.trg_ts_state            = cell2mat(tmp_trg_ts_state);

%%% For all coherence levels %%%
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
    out.mecc_state(iCoh)  	= nanmedian(cellfun(@(x) nanmedian(x(end-nSample:end)), in.js_ecc(cIdx)));
    out.ecc_state{iCoh}   	= cellfun(@(x) nanmedian(x(end-nSample:end)), in.js_ecc(cIdx));
   
    % Joystick accuracy [state-wise]
    tmp_js_dir_state        = in.js_dir(cIdx);
    tmp_rdp_dir_state       = in.rdp_dir(cIdx);
    
    clear js_acc_state js_dev
    for iState = 1:length(tmp_js_dir_state)
        js_dev                  = rad2deg(circ_dist(deg2rad(tmp_js_dir_state{iState}(end-nSample:end)),deg2rad(tmp_rdp_dir_state(iState))));  % Minimum RDP-Joystick difference
        js_acc_state(iState)	= nanmean(abs(1 - abs(js_dev) / 180));         	% Joystick accuracy
    end
    out.macc_state(iCoh)  	= nanmedian(js_acc_state);
    out.acc_state{iCoh}    	= js_acc_state';
    
    if length(out.acc_state{iCoh}) ~= length(out.ecc_state{iCoh})
        error('stop');
    end
    
    % Joystick accuracy [before first target]
    t1_outc                 = cellfun(@(x) x(1), in.trg_hit,'UniformOutput', false);
    t1_ts                   = cellfun(@(x) x(1), cellfun(@double,in.trg_ts,'UniformOutput', false));
    f1_ts                   = cellfun(@(x) x(1), cellfun(@double,in.frme_ts,'UniformOutput', false));
    trgIdx                  = (t1_ts-f1_ts) >= 1e6; % First target appearance 1s afer direction change
    rdp_dir                 = in.rdp_dir(cIdx & in.trg_shown & trgIdx);
    js_dir                  = in.js_dir(cIdx & in.trg_shown & trgIdx);
    js_ecc_tmp              = in.js_ecc(cIdx & in.trg_shown & trgIdx);
    frmes                   = in.frme_ts(cIdx & in.trg_shown & trgIdx);
    trg1_ts                 = t1_ts(cIdx & in.trg_shown & trgIdx);
    
    clear js_acc js_ecc js_ecc_end_of_state
    for iTrg = 1:length(rdp_dir)
        clear js_dev
        smpl_idx            = find(frmes{iTrg} < trg1_ts(iTrg),1,'last')-nSample : find(frmes{iTrg} < trg1_ts(iTrg),1,'last');
        js_dev              = rad2deg(circ_dist(deg2rad(js_dir{iTrg}(smpl_idx)),deg2rad(rdp_dir(iTrg))));  % Minimum RDP-Joystick difference
        js_acc(iTrg)        = nanmean(abs(1 - abs(js_dev) / 180));         	% Joystick accuracy
        js_ecc(iTrg)      	= nanmean(js_ecc_tmp{iTrg}(smpl_idx));
        js_ecc_end_of_state(iTrg) = nanmean(js_ecc_tmp{iTrg}(end-nSample:end));%%%
    end
    
    out.macc_trg(iCoh)      = nanmedian(js_acc);
    out.mecc_trg(iCoh)      = nanmedian(js_ecc);
    out.acc_trg{iCoh}     	= js_acc;
    out.ecc_trg{iCoh}       = js_ecc;
    out.ecc_end_of_state{iCoh} = js_ecc_end_of_state;%%%
    out.outc_trg{iCoh}      = cell2mat(t1_outc(cIdx & in.trg_shown & trgIdx));
    out.carr(iCoh)         	= snr(iCoh);
end
end

function plot_trg_response(in)

outc = logical(cell2mat(in.trg_hit(logical(in.trg_shown))'));
ecc = cell2mat(in.trg_ecc(logical(in.trg_shown))');
acc = cell2mat(in.trg_acc(logical(in.trg_shown))');
acc_indx = acc > .9;
edges = 0:.05:1;

figure
subplot(2,2,1)
histogram(acc(outc),edges)
hold on
histogram(acc(~outc),edges)
xlabel('Accuracy')
ylabel('# Targets')
title('Accuracy')
legend('Hit', 'Miss', 'Location', 'northwest')

subplot(2,2,2)
histogram(ecc(outc),edges)
hold on
histogram(ecc(~outc),edges)
xlabel('Eccentricity')
ylabel('# Targets')
title('Eccentricity')

subplot(2,2,3)
histogram(ecc(outc & acc_indx),edges)
hold on
histogram(ecc(outc & ~acc_indx),edges)
xlabel('Eccentricity')
ylabel('# Targets')
title('Eccentricity [Hits]')
legend('Acc > .9', 'Acc < .9', 'Location', 'northwest')

subplot(2,2,4)
histogram(ecc(~outc & acc_indx),edges)
hold on
histogram(ecc(~outc & ~acc_indx),edges)
xlabel('Eccentricity')
ylabel('# Targets')
title('Eccentricity [Misses]')
end
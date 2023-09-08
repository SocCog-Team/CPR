clear all
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/

% Import subject summary spreadsheet
pth                    	= '/Volumes/T7_Shield/CPR_psychophysics/';      % Local hard drive
dest_dir            	= '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/Stats/';
excel_sheet             = readtable([pth 'Subjects_summary.xlsx']);     % Spreadsheet
sbj_lst               	= excel_sheet.Abbreviation;                   	% Subject ID list
nSample             	= 29;                                           % Last 30 sample [end-nSample]
tbl                     = [];                                           % Initialise

% For each subject...
for iSubj = 1:length(sbj_lst)
    disp(['Processing: ' sbj_lst{iSubj}])
    
    % Go to relevant directory
    raw_pth          	= [pth sbj_lst{iSubj} '/raw/'];            
    cd(raw_pth)
    mwk2_files         	= dir('*.mwk2');

    data_pth          	= [pth sbj_lst{iSubj} '/summary/'];
    cd(data_pth)  
    mat_files         	= dir('*.mat');                                 
    
    % For each file...
    for iFile = 1:length(mat_files) 
        if contains(mat_files(iFile).name,'CPRsolo') || contains(mat_files(iFile).name,'CPRagent')            
            
            % Load pre-processed data table
            tmp_solo   	= load(mat_files(iFile).name); 
            
            % Exclude sessions with different coherence level
            if sum(unique(tmp_solo.t.rdp_coh) < .1) > 2
                continue
            end
            
            % Get recording time
            clear fIdx
            for iRaw = 1:length(mwk2_files)
                tmp       	= split(mat_files(iFile).name,'_');
                
                if contains(mat_files(iFile).name,'CPRagent')
                    tmp{2}  = lower(sbj_lst{iSubj});
                end
                
                tmp         = join(tmp(1:4),'_');
                fIdx(iRaw)  = contains(mwk2_files(iRaw).name, tmp);
            end
 
            % Build and concatenate data table
            out       	= build_stats_table(tmp_solo.t,excel_sheet,nSample,lower(sbj_lst{iSubj}),mwk2_files(fIdx).date);
            tbl         = [tbl;out];
        end
    end
end

% For each dyad....
for iDyad = 19:71
    disp(['Processing Dyad: ' num2str(iDyad)])

    % Exclude sessions
    if iDyad == 24 || iDyad == 36 || iDyad == 55
        continue
    end
    
    % Go to relevant directory
    cd([pth ['Dyad' num2str(iDyad)] '/raw/'])
    mwk2_files         	= dir('*.mwk2');
    
    % Go to relevant directory
    cd([pth ['Dyad' num2str(iDyad)] '/summary/'])
    mat_files              	= dir('*.mat');

    % Extract subject ID
    tmp_fid                 = split(mwk2_files(1).name,'_');
    id_dyad                 = [{tmp_fid{2}(1:3)}, {tmp_fid{2}(4:end)}];
    
    % For each file...
    for iFile = 1:length(mat_files)
        
        % Determine recording date
        clear fIdx
        for iRaw = 1:length(mwk2_files)
            tmp       	= split(mat_files(iFile).name,'_');
            fIdx(iRaw)  = contains(mwk2_files(iRaw).name, tmp{4});
        end  
        
        % Load pre-processed data table
        tmp_dyad   	= load(mat_files(iFile).name);  
        
        % Determine second player
        other_plyer = id_dyad{~(strcmp(tmp_dyad.t.ID(1), id_dyad))};
            
        % Build and concatenate data table
        out       	= build_stats_table(tmp_dyad.t,excel_sheet,nSample,other_plyer,mwk2_files(fIdx).date);
        tbl         = [tbl;out];
    end
end

save([dest_dir 'Stats_tbl.mat'], 'tbl', '-v7.3');
writetable(tbl,[dest_dir 'Stats_tbl.csv']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = build_stats_table(in,excel_sheet,nSample,human,dte)

% Experiment + subject parameter
tbl_0                  = in(:,1);

if strcmp(in.exp(1),'CPRdyadic') || (strcmp(in.exp(1),'CPRagent') && strcmp(in.ID(1),'agnt'))
    partner          	= repmat(convertCharsToStrings(lower(human)),[size(tbl_0,1),1]);
elseif strcmp(in.exp(1),'CPRagent') && strcmp(in.ID(1),human)
    partner            	= repmat(convertCharsToStrings('computer'),[size(tbl_0,1),1]);
elseif strcmp(in.exp(1),'CPRsolo')
    partner            	= repmat(convertCharsToStrings('none'),[size(tbl_0,1),1]);    
end

if strcmp(tbl_0.ID(1),'agnt')
    gender              = repmat({'none'},[size(tbl_0,1),1]);
    age              	= repmat(999,[size(tbl_0,1),1]);
    hand               	= repmat({'none'},[size(tbl_0,1),1]);
else
    sIdx              	= cellfun(@(x) strcmp(lower(x),tbl_0.ID(1)),excel_sheet.Abbreviation);
    gender           	= repmat(excel_sheet.Gender(sIdx),[size(tbl_0,1),1]);
    age                	= repmat(excel_sheet.Age(sIdx),[size(tbl_0,1),1]);
    hand               	= repmat(excel_sheet.Handedness(sIdx),[size(tbl_0,1),1]);    
end

ts_rec                  = repmat(datetime(dte),[size(tbl_0,1),1]);
tbl_1                   = [table(ts_rec,'VariableNames',{'datetime'}) in(:,3) table(partner,'VariableNames',{'other_player'}) in(:,4:6) in(:,8:12)];
tbl_2                   = table(gender,age,hand,'VariableNames',{'gender','age','hand'});

% Number of targets
ntrg                    = cellfun(@length,in.trg_hit);
ntrg(in.trg_shown == 0) = 0;

% Number of hits
nhits                   = cellfun(@sum,in.trg_hit);

% Target score
tbl_3                    = in(:,17);

% Remove NaNs
for iState = 1:height(tbl_3)
    for iTrg = 1:length(tbl_3{iState,1}{1})
        nanIdx          = isnan(tbl_3{iState,1}{1});
        if sum(nanIdx)>0
            tbl_3{iState,1}{1}(nanIdx) = 0;
        end
    end
end

% Joystick displacement
try
    js_ecc          	= cellfun(@(x) nanmedian(x(end-nSample:end)), in.js_ecc);
catch
    for iState = 1:length(in.js_ecc)
        if length(in.js_ecc{iState}) > 100
            js_ecc(iState,1)  	= nanmedian(in.js_ecc{iState}(end-nSample:end));
        else
            js_ecc(iState,1)	= nan;
        end
    end  
end

% Joystick accuracy [before first target]
t1_ts                   = cellfun(@(x) x(1), in.trg_ts);
f1_ts                   = cellfun(@(x) x(1), in.frme_ts);
trgIdx                  = (t1_ts-f1_ts) >= 1e6; % 1sec

for iState = 1:length(trgIdx)
    if trgIdx(iState)
        clear frmes js_dev
        
        frmes         	= in.frme_ts{iState};
        smpl_idx      	= find(frmes < t1_ts(iState),1,'last')-nSample : find(frmes < t1_ts(iState),1,'last');
        
        rdp_dir       	= deg2rad(in.rdp_dir(iState));
        js_dir        	= deg2rad(in.js_dir{iState}(smpl_idx));
        
        js_dev        	= rad2deg(circ_dist(js_dir,rdp_dir));           % Minimum RDP-Joystick difference
        js_acc(iState)	= nanmean(abs(1 - abs(js_dev) / 180));         	% Joystick accuracy
    else
        js_acc(iState)	= nan;
    end
end

tbl_4                   = table(ntrg,nhits,'VariableNames',{'n_targets','n_hits'});
tbl_5                   = table(js_ecc,js_acc', 'VariableNames',{'js_eccentricity', 'js_accuracy'});

out                     = [tbl_0 tbl_2 tbl_1 tbl_4 tbl_3 tbl_5]; % Output table
end
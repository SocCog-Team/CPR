function t = CPR_construct_table(subj, fname, d, trl, idx, write_file)

% This function constructs a state-wise data table based on a previously 
% imported .mwk2 data structure.
%
% Input:  	.subj           String, Subject ID
%           .fid            String, File identifier
%           .d              Structure, Contains raw data
%          	.trl            Structure, Contains trial information
%           .idx            Structure, Contains variable-wise indices
%           .write_file     Boolean, Indicates if table is to be saved
%
% Output:   .t              Table, State-wise data table containing 
%                           stimulus parameters as well as behavioural 
%                           responses 
%
% Example:	t = CPR_construct_table('fxs','20210504',d, trl, idx, true)
%
% Known bugs:
%
% Feature wish list:
%
% Version history
%   1.0     (fxs 2021-05-04) Initial version.
%   1.1     (fxs 2021-05-25) Added compatibility with dyadic setting.
%   1.2     (fxs 2021-11-24) Compatibility with new MWEL paradigm.
%   1.2     (fxs 2021-11-24) Frame-wise output.                            

%% Initialise

% Table specs
cc    	= 0;
var     = {'ID', 'string'; ...              % Subject identity
    'date', 'string'; ...                   % Experiment date
    'exp', 'string'; ...                    % Trial/Block number
    'setup', 'string'; ...                  % Setup
    'block', 'double'; ...                  % Experimental block
    'trl_no', 'double'; ...                 % Trial/Block number
    'trl_dur', 'double'; ...                % Trial/Block duration
    'ss_no', 'double'; ...                  % Steady state number
    'ss_coh', 'double'; ...                 % Steady state coherence
    'ss_dur', 'double'; ...                 % Steady state duration
    'rdp_dir', 'double';...                 % RDP direction
    'trg_shown', 'double'; ...              % Target shown?
    'trg_hit', 'cell'; ...                  % Target hit?
    'trg_score', 'cell'; ...                % Target score
    'trg_ts', 'cell'; ...                   % Target timestamp
    'frme_ts', 'cell'; ...                  % Frame timestamp
    'js_dir', 'cell'; ...                   % Joystick direction
    'js_str', 'cell'; ...                   % Joystick strength
    'trl_frme_ts', 'cell';...               % Trial-wise frame timestamp
    'trl_rdp_dir', 'cell'; ...              % Trial-wise RDP direction
    'trl_rdp_coh', 'cell'; ...              % Trial-wise RDP coherence
    'trl_js_dir', 'cell'; ...               % Trial-wise joystick direction
    'trl_js_str', 'cell'; ...               % Trial-wise joystick strength
    'trl_eye_x', 'cell'; ...                % Trial-wise eye data: x
    'trl_eye_y', 'cell'; ...                % Trial-wise eye data: y
    'trl_eye_ts', 'cell';...                % Trial-wise eye data: timestamps
    'trl_fix', 'cell';...                   % Trial-wise fixation breaks
    'trl_fix_ts', 'cell'};                  % Trial-wise fixation break timestamps

% Initialise table
t       = table('Size',[5000, size(var,1)],...
    'VariableTypes',var(:,2),...
    'VariableNames',var(:,1));

%% Build table

fprintf('Building table...\n')

% Extract title information
fid = strsplit(fname,'_');
psy4 = strcmp(fid{2}(1:3),subj);    

if psy4 == 1
    setup = 'psycho4';
else
    setup = 'psycho3';
end

% Trial loop
for iTrl = 1:length(trl.tEnd)
    
    trlIdx                  = [];
    trlIdx                  = d.time >= trl.tOn(iTrl) & d.time <= trl.tEnd(iTrl);   	% Trial index
    trl.state_ts{iTrl}      = getTrialData(d.time, trlIdx, idx.steady_duration);        % State ON duration timestamps
    trl.ssdur{iTrl}         = getTrialData(d.value, trlIdx, idx.steady_duration);       % Steady state duration for this trials
    trl.coh{iTrl}           = getTrialData(d.value, trlIdx, idx.RDP_coh);               % Trial RDP coherence level
    trl.dir{iTrl}           = getTrialData(d.value, trlIdx, idx.RDP_dir);               % Trial RDP direction
    
    if length(trl.ssdur{iTrl}) < 2
        continue
    end
    
    % Steady state loop
    for iSS = 1:length(trl.ssdur{iTrl})
        
        % Steady state index
        if iSS < length(trl.ssdur{iTrl})
            ssIdx           = d.time >= trl.state_ts{iTrl}(iSS) & d.time < trl.state_ts{iTrl}(iSS+1); 
        else
            ssIdx           = d.time >= trl.state_ts{iTrl}(iSS) & d.time <= trl.tEnd(iTrl);
        end
        
        % Fill in table: Trial/State/Stimulus parameter
        cc                  = cc+1;                                                     % Steady state counter
        t.ID{cc}            = subj;                                                     % Subject ID
        t.date(cc)          = fid{1};                                                   % Experimental date
        t.exp{cc}           = fid{3};                                                   % Experimental condition
        t.setup{cc}         = setup;                                                    % Setup
        t.block(cc)         = str2double(fid{4}(end));                                  % Experimental block
        t.trl_no(cc)        = iTrl;                                                     % Trial number
        t.trl_dur(cc)       = (trl.tEnd(iTrl) - trl.tOn(iTrl)) ./ 1e3;                  % Trial duration [ms]
        t.ss_no(cc)         = cc;                                                       % Steady state counter  
        t.ss_coh(cc)        = trl.coh{iTrl}(iSS);                                       % Steady state coherence
        t.rdp_dir(cc)       = mod(trl.dir{iTrl}(iSS),360);                              % Stimulus direction
        t.ss_dur(cc)        = getTrialData(d.value, ssIdx, idx.steady_duration);      	% Steady state duration
        t.trg_shown(cc)     = iscell(getTrialData(d.value, ssIdx, idx.outcome));        % Target shown?
        t.trg_hit{cc}       = strcmp(getTrialData(d.value, ssIdx, idx.outcome), 'hit'); % Target collected?
        t.trg_score{cc}     = nan(1,length(t.trg_hit{cc}));                             % Target score vector
        tmp_trg_val        	= getTrialData(d.value, ssIdx, idx.trg);                    % Target trigger
        tmp_trg_ts        	= getTrialData(d.time, ssIdx, idx.trg);                     % Target timestamp
        
        if sum(tmp_trg_val) > 0
            t.trg_ts{cc}                        = tmp_trg_ts(logical(tmp_trg_val));
            if sum(t.trg_hit{cc}) > 0
                tmp_trg_score                   = getTrialData(d.value, ssIdx, idx.reward);           % Target score
                t.trg_score{cc}(t.trg_hit{cc})  = tmp_trg_score;
            end
        else
            t.trg_ts{cc}                        = nan;
        end
        
        % Fill in table: Behavioural response
        % (1) For each steady state...
        
        % Tmp variables based on joystick sampling rate
        tmp_js_ts{cc}       = getTrialData(d.time, ssIdx, idx.JS_str);                  % Timestamps: Joystick strength
        tmp_js_dir{cc}      = getTrialData(d.value, ssIdx, idx.JS_dir);                 % Joystick direction
        tmp_js_str{cc}      = getTrialData(d.value, ssIdx, idx.JS_str);                 % Joystick strength
        
        t.frme_ts{cc}       = getTrialData(d.time, ssIdx, idx.frame);                   % Frame timestamps
        
        % Build frame-wise vector for joystick data
        for iFrme = 1:length(t.frme_ts{cc})
            fIdx                        = [];
            fIdx                        = find(tmp_js_ts{cc} < t.frme_ts{cc}(iFrme),1,'last');  % Extract last entry before frame onset  
            if sum(fIdx) == 0 || isempty(fIdx) || fIdx > length(tmp_js_dir{cc}) || fIdx > length(tmp_js_str{cc})        
                t.js_dir{cc}(iFrme)     = nan;                                          % Write to vector
                t.js_str{cc}(iFrme)     = nan;
            else
                t.js_dir{cc}(iFrme)     = tmp_js_dir{cc}(fIdx);                           
                t.js_str{cc}(iFrme)     = tmp_js_str{cc}(fIdx);
            end
        end
             
        % (2) For entire trial...
        if iSS == 1
            
            % TMP trial variables
            tmp_trl_js_ts{cc}     	= getTrialData(d.time, trlIdx, idx.JS_str);
            tmp_trl_js_dir{cc}   	= mod(getTrialData(d.value, trlIdx, idx.JS_dir),360);
            tmp_trl_js_str{cc}    	= getTrialData(d.value, trlIdx, idx.JS_str);
          	tmp_trl_rdp_dir_ts{cc}  = getTrialData(d.time, trlIdx, idx.RDP_dir);
            tmp_trl_rdp_dir{cc}   	= mod(getTrialData(d.value, trlIdx, idx.RDP_dir),360);
            tmp_trl_rdp_coh_ts{cc}	= getTrialData(d.time, trlIdx, idx.RDP_coh);
            tmp_trl_rdp_coh{cc}    	= getTrialData(d.value, trlIdx, idx.RDP_coh);
     
            % Frame timestamps
            t.trl_frme_ts{cc}       =  getTrialData(d.time, trlIdx, idx.frame);       	

            % Frame-wise RDP direction
            for iFrme = 1:length(t.trl_frme_ts{cc})
                fdIdx                           = [];
                fdIdx                           = find(tmp_trl_rdp_dir_ts{cc} < t.trl_frme_ts{cc}(iFrme),1,'last'); 
                
                if sum(fdIdx) == 0 || isempty(fdIdx) || fdIdx > length(tmp_trl_rdp_dir{cc})
                    t.trl_rdp_dir{cc}(iFrme) 	= nan;                           	
                else
                    t.trl_rdp_dir{cc}(iFrme)  	= tmp_trl_rdp_dir{cc}(fdIdx);                           
                end
            end
            
            % Frame-wise RDP coherence
            for iFrme = 1:length(t.trl_frme_ts{cc})
                fcIdx                           = [];
                fcIdx                           = find(tmp_trl_rdp_coh_ts{cc} < t.trl_frme_ts{cc}(iFrme),1,'last');
                
                if sum(fcIdx) == 0 || isempty(fcIdx) || fcIdx > length(tmp_trl_rdp_coh{cc})
                    t.trl_rdp_coh{cc}(iFrme)	= nan;
                else
                    t.trl_rdp_coh{cc}(iFrme)	= tmp_trl_rdp_coh{cc}(fcIdx);
                end
            end
            
            % Frame-wise joystick data
            for iFrme = 1:length(t.trl_frme_ts{cc})
                fjIdx                       = [];
                fjIdx                       = find(tmp_trl_js_ts{cc} < t.trl_frme_ts{cc}(iFrme),1,'last');  
                
                if sum(fjIdx) == 0 || isempty(fjIdx) || fjIdx > length(tmp_trl_js_dir{cc}) || fjIdx > length(tmp_trl_js_str{cc})  
                    t.trl_js_dir{cc}(iFrme)     = nan;                           	% Write to vector
                    t.trl_js_str{cc}(iFrme)     = nan;
                else
                    t.trl_js_dir{cc}(iFrme)     = tmp_trl_js_dir{cc}(fjIdx);                           	% Write to vector
                    t.trl_js_str{cc}(iFrme)     = tmp_trl_js_str{cc}(fjIdx);
                end
            end
            
            t.trl_eye_x{cc}         = getTrialData(d.value, trlIdx, idx.eye_x_dva);
            t.trl_eye_y{cc}         = getTrialData(d.value, trlIdx, idx.eye_y_dva);
            t.trl_eye_ts{cc}        = getTrialData(d.time, trlIdx, idx.eye_y_dva);
            t.trl_fix{cc}           = getTrialData(d.value, trlIdx, idx.fixation);
            t.trl_fix_ts{cc}        = getTrialData(d.time, trlIdx, idx.fixation);
        else
            t.trl_frme_ts{cc}       = nan;
            t.trl_rdp_dir{cc}       = nan;
            t.trl_rdp_coh{cc}       = nan;
            t.trl_js_dir{cc}        = nan;
            t.trl_js_str{cc}        = nan;
            t.trl_eye_x{cc}         = nan;
            t.trl_eye_y{cc}         = nan;
            t.trl_eye_ts{cc}        = nan;
            t.trl_fix{cc}           = nan;
            t.trl_fix_ts{cc}        = nan;
        end
    end
end

t(ismissing(t.ID),:)                = [];
fprintf('Done!\n')

if write_file
    fprintf(['Save table for subject ' subj '...\n'])
    pth_raw = pwd;
    cd ../summary
    
    % Save as .mat file
    save([fid{1} '_' subj '_' fid{3} '_' fid{4} '_tbl.mat'], 't', '-v7.3')
    
    cd(pth_raw)
    fprintf('Done!\n')
end
end
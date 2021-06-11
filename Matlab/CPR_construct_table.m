function t = CPR_construct_table(subj, fid, d, trl, idx, write_file)

% This function constructs a state-wise data table based on a previously 
% imported .mwk2 data structure.
%
% Input:  	.subj           String, Subject ID
%           .fid            String, File identifier
%           .d              Structure, Contains raw data
%          	.trl            Structure, Contains trial information
%           .idx            Structure, Contains variable-wise indices
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
%   1.0     (fxs 2020-05-04) Initial version.

%% Initialise

% Table specs
cc    	= 0;
var     = {'ID', 'string'; ...
    'trl_no', 'double'; ...
    'trl_dur', 'double'; ...
    'ss_no', 'double'; ...
    'ss_coh', 'double'; ...
    'ss_dur', 'double'; ...
    'rdp_dir', 'double';...
    'trg_shown', 'logical'; ...
    'trg_hit', 'logical'; ...
    'trg_ts', 'double'; ...
    'js_dir', 'cell'; ...
    'js_str', 'cell'; ...
    'js_ts', 'cell'; ...   
    'trl_rdp_dir', 'cell'; ...
    'trl_rdp_dir_ts', 'cell'; ...
    'trl_rdp_coh', 'cell'; ...
    'trl_rdp_coh_ts', 'cell'; ...
    'trl_js_dir', 'cell'; ...
    'trl_js_str', 'cell'; ...
    'trl_js_ts', 'cell';...
    'trl_eye_x', 'cell'; ...
    'trl_eye_y', 'cell'; ...
    'trl_eye_ts', 'cell';...
    'trl_fix', 'cell';...
    'trl_fix_ts', 'cell'};


% Initialise table
t       = table('Size',[5000, size(var,1)],...
    'VariableTypes',var(:,2),...
    'VariableNames',var(:,1));

%% Build table

% Trial loop
for iTrl = 1:length(trl.tEnd)
    
    trlIdx                  = [];
    trlIdx                  = d.time >= trl.tOn(iTrl) & d.time <= trl.tEnd(iTrl);   	% Trial index
    trl.dir_ts{iTrl}        = getTrialData(d.time, trlIdx, idx.stateOn);             	% Timestamps of RDP direction changes [Steady state start, Ignore first value -> randomly drawn direction]
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
            ssIdx           = d.time >= trl.dir_ts{iTrl}(iSS+1) & d.time < trl.dir_ts{iTrl}(iSS+2); % 1st timestamp irrelevant
        else
            ssIdx           = d.time >= trl.dir_ts{iTrl}(iSS+1) & d.time <= trl.tEnd(iTrl);
        end
        
        % Fill in table: Trial/State/Stimulus parameter
        cc                  = cc+1;                                                     % Steady state counter
        t.ID{cc}            = [subj '_' fid];                                           % Subject ID
        t.trl_no(cc)        = iTrl;                                                     % Trial number
        t.trl_dur(cc)       = (trl.tEnd(iTrl) - trl.tOn(iTrl)) ./ 1e3;                  % Trial duration [ms]
        t.ss_no(cc)         = cc;                                                       % Steady state counter
        
        if strcmp(subj, 'cla') || strcmp(subj, 'nil') || str2num(fid) < 20210401
            t.ss_coh(cc)  	= trl.coh{iTrl}(end);
            t.rdp_dir(cc)  	= mod(trl.dir{iTrl}(iSS),360);
        else
            t.ss_coh(cc)   	= getTrialData(d.value, ssIdx, idx.RDP_coh);                % Steady state coherence
            t.rdp_dir(cc)   = mod(getTrialData(d.value, ssIdx, idx.RDP_dir),360);       % Stimulus direction
        end
        
        t.ss_dur(cc)        = trl.ssdur{iTrl}(iSS);                                     % Steady state duration
        t.trg_shown(cc)     = iscell(getTrialData(d.value, ssIdx, idx.outcome));        % Target shown?
        t.trg_hit(cc)       = strcmp(getTrialData(d.value, ssIdx, idx.outcome), 'hit'); % Target collected?
        trg_val             = getTrialData(d.value, ssIdx, idx.trg);                    % Target trigger
        trg_ts              = getTrialData(d.time, ssIdx, idx.trg);                     % Target timestamp
        
        if sum(trg_val) == 1
            t.trg_ts(cc)  	= trg_ts(logical(trg_val));
        else
            t.trg_ts(cc) 	= nan;
        end
        clear trg_val trg_ts
        
        % Fill in table: Behavioural response
        % (1) For each steady state...
        t.js_dir{cc}        = getTrialData(d.value, ssIdx, idx.JS_dir);                 % Joystick direction
        t.js_dir_ts{cc}     = getTrialData(d.time, ssIdx, idx.JS_dir);                  % Timestamps: Joystick direction
        t.js_str{cc}        = getTrialData(d.value, ssIdx, idx.JS_str);                 % Joystick strength
        t.js_str_ts{cc}     = getTrialData(d.time, ssIdx, idx.JS_str);                  % Timestamps: Joystick strength
        
        % (2) For entire trial...
        if iSS == 1
            t.trl_rdp_dir{cc}       = mod(getTrialData(d.value, trlIdx, idx.RDP_dir),360);
            t.trl_rdp_dir_ts{cc}    = getTrialData(d.time, trlIdx, idx.RDP_dir);
            t.trl_rdp_coh{cc}       = getTrialData(d.value, trlIdx, idx.RDP_coh);
            t.trl_rdp_coh_ts{cc}    = getTrialData(d.time, trlIdx, idx.RDP_coh);
            t.trl_js_dir{cc}        = mod(getTrialData(d.value, trlIdx, idx.JS_dir),360);
            t.trl_js_str{cc}        = getTrialData(d.value, trlIdx, idx.JS_str);
            t.trl_js_ts{cc}         = getTrialData(d.time, trlIdx, idx.JS_str);
            t.trl_eye_x{cc}         = getTrialData(d.value, trlIdx, idx.eye_x_dva);
            t.trl_eye_y{cc}         = getTrialData(d.value, trlIdx, idx.eye_y_dva);
            t.trl_eye_ts{cc}        = getTrialData(d.time, trlIdx, idx.eye_y_dva);
            t.trl_fix{cc}           = getTrialData(d.value, trlIdx, idx.fixation);
            t.trl_fix_ts{cc}        = getTrialData(d.time, trlIdx, idx.fixation);
        else
            t.trl_rdp_dir{cc}       = nan;
            t.trl_rdp_dir_ts{cc}    = nan;
            t.trl_rdp_coh{cc}       = nan;
            t.trl_rdp_coh_ts{cc}    = nan;
            t.trl_js_dir{cc}        = nan;
            t.trl_js_str{cc}        = nan;
            t.trl_js_ts{cc}         = nan;
            t.trl_eye_x{cc}         = nan;
            t.trl_eye_y{cc}         = nan;
            t.trl_eye_ts{cc}        = nan;
            t.trl_fix{cc}           = nan;
            t.trl_fix_ts{cc}        = nan;
        end
    end
    
    % Temporary fix: coherence per steady state
    if strcmp(subj, 'cla') || strcmp(subj, 'nil') || str2num(fid) < 20210401
    else
        b = repmat(trl.coh{iTrl},10,1);
        if iTrl == 1
            c = reshape(b(:,2:end),size(b,1)*(size(b,2)-1),1);
        else
            c = reshape(b,size(b,1)*(size(b,2)),1);
        end
        strt = (iTrl*100)-99;
        t.ss_coh(strt:strt+length(c)-1) = c;
    end
end

t(ismissing(t.ID),:)                = [];

if write_file
    disp('Save table...')
    save([subj '_' fid '_tbl.mat'], 't', '-v7.3')                                                   % Save as .mat file
    disp('Done!')
end
end
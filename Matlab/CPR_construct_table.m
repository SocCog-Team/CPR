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
%   1.3     (fxs 2022-09-05) Block-wise target processing.

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
    'ss_dur', 'double'; ...                 % Steady state duration
    'rdp_coh', 'double'; ...                % RDP coherence
    'rdp_dir', 'double';...                 % RDP direction
    'trg_shown', 'double'; ...              % Target shown?
    'trg_ts', 'cell'; ...                   % Target timestamp
    'trg_hit', 'cell'; ...                  % Target hit?
    'trg_acc', 'cell'; ...                  % Target accuracy
    'trg_ecc', 'cell'; ...                  % Target eccentricity
    'trg_score', 'cell'; ...                % Target score
    'frme_ts', 'cell'; ...                  % Frame timestamp
    'js_dir', 'cell'; ...                   % Joystick direction
    'js_ecc', 'cell'; ...                   % Joystick direction
    'rdp_dot', 'cell'};                     % RDP dot position

% Initialise table
t       = table('Size',[5000, size(var,1)],...
    'VariableTypes',var(:,2),...
    'VariableNames',var(:,1));

%% Build table

fprintf('Building table...\n')

% Extract title information
fid          	= strsplit(fname,'_');

if strcmp(fid{5},'psycho4') || strcmp(fid{5},'psycho3') || strcmp(fid{5},'physio4') || strcmp(fid{5},'physio3')
    setup   	= fid{5};
else
    psy4 = strcmp(fid{2}(1:3),subj);
    
    if psy4 == 1
        setup 	= 'psycho4';
    else
        setup 	= 'psycho3';
    end
end

% Stimulus cycle (i.e. trial) loop
for iTrl = 1:length(trl.tEnd)
    
    % Trial index
    trlIdx                  = [];
    trlIdx                  = d.time >= trl.tOn(iTrl) & d.time <= trl.tEnd(iTrl);
    
    % Extract dots
    tmp_dp                	= d.value(idx.RDP_dot);
    tmp_dp_ts           	= d.time(idx.RDP_dot);
    tmp_frme_ts             = d.time(idx.frame);
    
    % Stimulus trial data
    tmp_dir_ts              = getTrialData(d.time, trlIdx, idx.RDP_dir);                    % RDP_direction timestamps
    tmp_dir                 = getTrialData(d.value, trlIdx, idx.RDP_dir);                   % RDP_direction
    tmp_coh                 = getTrialData(d.value, trlIdx, idx.RDP_coh);                   % RDP coherence
    tmp_coh_ts            	= getTrialData(d.time, trlIdx, idx.RDP_coh);                    % RDP coherence timestamps
    
    % Remove first entries -> similar to last stimulus state at RDP_onset
    trl.state_ts{iTrl}      = tmp_dir_ts(2:end);                                            % State ON duration timestamps
    trl.state_dur{iTrl}    	= diff([tmp_dir_ts(2:end) trl.tEnd(iTrl)]) ./ 1e3;              % Steady state duration for this trials [ms]
    trl.dir{iTrl}           = tmp_dir(2:end);                                               % Trial RDP direction
    
    % Stimulus state loop
    for iSS = 1:length(trl.state_dur{iTrl})
        
        % Steady state index
        if iSS < length(trl.state_dur{iTrl})
            ssIdx           = d.time >= trl.state_ts{iTrl}(iSS) & d.time < trl.state_ts{iTrl}(iSS+1);
        else
            ssIdx           = d.time >= trl.state_ts{iTrl}(iSS) & d.time <= trl.tEnd(iTrl);
        end
        
        % Trial/State/Stimulus parameter
        cc                  = cc+1;                                                         % Steady state counter
        t.ID{cc}            = subj;                                                         % Subject ID
        t.date(cc)          = fid{1};                                                       % Experimental date
        t.exp{cc}           = fid{3};                                                       % Experimental condition
        t.setup{cc}         = setup;                                                        % Setup
        t.block(cc)         = str2double(fid{4}(end));                                      % Experimental block
        t.trl_no(cc)        = iTrl;                                                         % Trial number
        t.trl_dur(cc)       = (trl.tEnd(iTrl) - trl.tOn(iTrl)) ./ 1e3;                      % Trial duration [ms]
        t.ss_no(cc)         = cc;                                                           % Stimulus state counter
        t.ss_dur(cc)        = trl.state_dur{iTrl}(iSS);                                     % Steady state duration
        t.rdp_coh(cc)       = tmp_coh(find(tmp_coh_ts <= trl.state_ts{iTrl}(iSS),1,'last'));% Stimulus state coherence
        t.rdp_dir(cc)       = mod(trl.dir{iTrl}(iSS),360);                                  % Stimulus direction
        
        % Behavioural response: TMP_ variables based on joystick sampling rate
        tmp_js_ts{cc}       = getTrialData(d.time, ssIdx, idx.JS_str);                      % Timestamps: Joystick strength
        tmp_js_dir{cc}      = getTrialData(d.value, ssIdx, idx.JS_dir);                     % Joystick direction
        tmp_js_ecc{cc}      = getTrialData(d.value, ssIdx, idx.JS_str);                     % Joystick strength
        
        t.frme_ts{cc}       = getTrialData(d.time, ssIdx, idx.frame);                       % Frame timestamps
        
        % Build frame-wise vector for joystick data
        for iFrme = 1:length(t.frme_ts{cc})
            fIdx                  	= [];
            fIdx                 	= find(tmp_js_ts{cc} < t.frme_ts{cc}(iFrme),1,'last');  % Extract last entry before frame onset
            if sum(fIdx) == 0 || isempty(fIdx) || fIdx > length(tmp_js_dir{cc}) || fIdx > length(tmp_js_ecc{cc})
                t.js_dir{cc}(iFrme)	= nan;                                                  % Write to vector
                t.js_ecc{cc}(iFrme)	= nan;
            else
                t.js_dir{cc}(iFrme)	= tmp_js_dir{cc}(fIdx);
                t.js_ecc{cc}(iFrme)	= tmp_js_ecc{cc}(fIdx);
            end
        end
        
        % Add dot positions
        t.rdp_dot{cc}               = d.value(ssIdx & idx.RDP_dot);
    end
end

% Crop table
t(ismissing(t.ID),:)                = [];

% Extract target presentation times
trg_ts                              = d.time(idx.trg_on);
trg_ts                              = trg_ts(logical(cell2mat(d.value(idx.trg_on))));

% Extract target outcome
outcome                             = d.value(idx.outcome);
outcome_ts                          = d.time(idx.outcome);
rew_score                           = d.value(idx.reward);                                      % Reward scores
rew_ts                              = d.time(idx.reward);                                       % Reward scores

%bExtract state and joystick parameter
state_on                            = cellfun(@(x) x(1), t.frme_ts);                            % State onset timestamps
js_dir_val                       	= d.value(idx.JS_dir);
js_ecc_val                       	= d.value(idx.JS_str);

% o = d.value(idx.outcome);
% tm = d.time(idx.outcome);
% for i = 1:length(trg_ts)
%     s(i) = o(find(trg_ts(i) < tm,1,'first'));
% end

%xxx = d.time(idx.reward)
% trg_ts(21:23)-outcome_ts(21:23)
% trg_ts(21:23)-rew_ts(21:23)

% Assign target data to stimulus states
for iTarget = 1:length(trg_ts)
    iState                       	= find(trg_ts(iTarget) >= state_on, 1,'last');              % Find target state
    
    if isempty(iState)
        continue
    end
    
    t.trg_shown(iState)             = true;                                                     % Indicate target state
    t.trg_ts{iState}                = [t.trg_ts{iState} trg_ts(iTarget)];                       % Extract target timestamp
    
    tmp_outcome                     = outcome(find(outcome_ts > trg_ts(iTarget),1,'first'));    % Extract target outcome
    t.trg_hit{iState}               = [t.trg_hit{iState} , strcmp(tmp_outcome, 'hit')];         % Indicate if hit or miss
    
    smple                           = find(d.time(idx.JS_dir) <= trg_ts(iTarget), 1, 'last');  	% Extract last JS sample before target
    js_dev                          = rad2deg(circ_dist(deg2rad(js_dir_val{smple}),deg2rad(t.rdp_dir(iState)))); % Get circular distance to RDP direction
    js_acc                       	= abs(1 - abs(js_dev / 180));                               % Calculate accuracy
    
    t.trg_acc{iState}               = [t.trg_acc{iState} js_acc];                               % Add accuracy to vector
    t.trg_ecc{iState}               = [t.trg_ecc{iState} js_ecc_val{smple}];                    % Add eccentricity to vector
    
     %%% THIS SEEMS PROBLEMATIC> >>> BUG!
     %%% Miss label but reward score ?!
     %%% ignore and take last score reading as fix?
    if strcmp(tmp_outcome, 'hit')
        rew_idx                     = find(rew_ts > trg_ts(iTarget),1,'first');                 % Get reward score
        last_score                  = rew_score{rew_idx-1};
        
        if iTarget == 1
            t.trg_score{iState}    	= [t.trg_score{iState} rew_score{rew_idx}-0];               % Add score to vector
        else
            t.trg_score{iState}   	= [t.trg_score{iState} rew_score{rew_idx}-rew_score{rew_idx-1}]; 
        end
        
    else
        t.trg_score{iState}       	= [t.trg_score{iState} 0];                                  % Add zero for missed targets
    end
end

% Fill in states without targets
no_trg_states                       = find(cellfun(@isempty,t.trg_ts));
for j = 1:length(no_trg_states)
    iState                          = no_trg_states(j);
    t.trg_shown(iState)          	= false;
    t.trg_ts{iState}              	= nan;
    t.trg_ecc{iState}              	= nan;
    t.trg_acc{iState}              	= nan;
    t.trg_hit{iState}             	= false;
    t.trg_score{iState}         	= nan;
end

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
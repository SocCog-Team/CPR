function [d,t] = PHY_daily_summary_test(data_pth, fname, cfg_pth, writeH5_flag, dest_dir_tbl, dest_dir_plot)
% Produces daily NHP training summary
%
% Input:  	.data_pth        	String, Data directory
%           .fname              String, File name
%                               Format: 'YYYMMDD_nil_CPRsolo_block1_physio4.mwk2'
%          	.cfg_pth            String, CFG file location
%           .writeH5_flag     	Boolean, True -> Write file
%           .dest_dir_tbl       String, Destination direcotry for data table
%           .dest_dir_plot      String, Destination direcotry for summary plot
%
% Output:   .t                  Table, State-wise data table containing
%                               stimulus parameters as well as behavioural
%                               responses
%           .d                  Imported MWorks data structure
%
% Example:	[d,t] = PHY_daily_summary('/Documents/MWorks/Data/', ...
%               '20240806_nil_CPRsolo_block1_physio4.mwk2', ...
%               'ownCloud/CPR_files/NHP_training/felix_nhp_solo.cfg', ...
%               false, 'ownCloud/CPR_files/NHP_training/mat/',...
%               'ownCloud/CPR_files/NHP_training/png/')
%
% Known bugs:
%
% Feature wish list:
%
% Version history
%   1.0     (fxs 2021-05-04) Initial version.

% Check input
if nargin < 6
    dest_dir_plot 	= '/Users/fschneider/ownCloud/CPR_files/NHP_training/png/';
end

if nargin < 5
    dest_dir_tbl    = '/Users/fschneider/ownCloud/CPR_files/NHP_training/mat/';
end

if nargin < 4
    writeH5_flag    = false;
end

if nargin < 3
    cfg_pth         = '/Users/fschneider/ownCloud/CPR_files/NHP_training/felix_nhp_solo.cfg';
end

if nargin < 2
    error('Specify data directory and file name!')
end

% Add relevant directories
addpath /Users/fschneider/ownCloud/CPR_files/NHP_training
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab/
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab

% Change into data directory
cd(data_pth)
close all

%% Check input arguments & import data

% Define variables for import
var_import  = {
    'INFO_', ...
    'TRIAL_', ...
    'IO_joystickDirection', ...
    'IO_joystickStrength',...
    'IO_fixation_flag',...
    'EYE_x_dva',...
    'EYE_y_dva',...
    '#stimDisplay'};
  
% Import (and H5-convert) data. If H5 exists already, set flag to false
d           = PHY_import_mwk2(fname, var_import, writeH5_flag, cfg_pth);

%% Construct data table

% Initialise variables
idx                         = [];
cyc                         = [];

% Create variable-specific indices
idx.cOn                     = d.event == 'TRIAL_start';
idx.cEnd                    = d.event == 'TRIAL_end';
idx.frame                   = d.event == 'STIM_displayUpdate';
idx.RDP_onset               = d.event == 'STIM_RDP_onset';
idx.RDP_dir                 = d.event == 'STIM_RDP_direction';
idx.RDP_coh                 = d.event == 'STIM_RDP_coherence';
idx.RDP_dot                	= d.event == 'STIM_RDP_dotPositions';
idx.trg_on                  = d.event == 'STIM_target_onset';
idx.JS_dir                  = d.event == 'IO_joystickDirection';
idx.JS_str                  = d.event == 'IO_joystickStrength';
idx.fixation              	= d.event == 'IO_fixation_flag';
idx.outcome                 = d.event == 'TRIAL_outcome';
idx.trg                     = d.event == 'TRIAL_reactionTrigger';
idx.eye_x_dva              	= d.event == 'EYE_x_dva';
idx.eye_y_dva             	= d.event == 'EYE_y_dva';
idx.reward                  = d.event == 'INFO_Juice_ml';
rew_str                     = 'ml';
    
% Get trial timestamps
cyc.cOn                     = d.time(idx.cOn);
cyc.cEnd                    = d.time(idx.cEnd);

% Extract title information
% FORMAT: fname = 'YYYMMDD_nil_CPRsolo_block1_physio4'

% if exist([dest_dir_tbl fname(1:end-5) '_tbl.mat'])
%     tmp = load([dest_dir_tbl fname(1:end-5) '_tbl.mat']);
%     t = tmp.t;
% else
    fid                         = strsplit(fname,'_');
    t                           = PHY_construct_table(fid{2}, fname, d, cyc, idx, false); % Don't write table through this fuction...
    
    % ...save data table as .mat file here
    save([dest_dir_tbl fname(1:end-5) '_tbl.mat'], 't', '-v7.3')
% end

%% Analyse behavior

%%% Performance [time window analysis]
nSample                 = 29;
nLag                    = 150;
solo_perf               = response_readout(t, nSample);

% Extract stimulus cycles for correlation analysis
sc_cnt                  = 0; % Reset cycle counter
cyc_boundary        	= [0; find(diff(t.cyc_no) ~= 0)];
[tmp_solo, sc_cnt]   	= extractCycles(cyc_boundary, t,[], sc_cnt);

%%% Correlation analysis
[solo_cr,~]            = PHY_correlation_analysis(tmp_solo, nLag, false);

%% Plot data
lw                      = 2;
fsz                     = 12;       
colmn                   = linspace(.1,.8,4);
height                  = [.77 .5 .23];
dim                     = [.15 .2];

% Initialise figure
f                       = figure('units','centimeters','position',[0 0 40 30]);

% PERFORMANCE: HIT RATE 
ax0                  	= axes('Position', [colmn(1) height(1) dim]); hold on
b                       = bar([solo_perf.hir;1-solo_perf.hir]','stacked');
b(1).FaceColor          = [.5 .5 .5]; % Hit
b(2).FaceColor          = [.1 .1 .1]; % Miss
lg                      = legend('Hit','Miss');
lg.Location             = 'northwest';
ax0.XTickLabel       	= {solo_perf.carr};
ax0.XLabel.String    	= 'Performance';
ax0.YLabel.String    	= 'Avg. correlation coefficient';
ax0.FontSize         	= fsz;

% PERFORMANCE: Fluid/Score 
ax1                  	= axes('Position', [colmn(2) height(1) dim]); hold on
[dat,cond]              = transformData(solo_perf.trg_score);
col                     = cool(length(solo_perf.carr));
vl                      = violinplot(dat,cond);
vl                      = improveViolin(vl,ax1,col,solo_perf.carr,fsz);
ax1.YLabel.String       = 'Juice/Target [ml]';

% ACCURACY PER COHERENCE --> state-aligned
ax2                  	= axes('Position', [colmn(3) height(1) dim]); hold on
[dat,cond]              = transformData(solo_perf.acc_state);
vl                      = violinplot(dat,cond);
vl                      = improveViolin(vl,ax2,col,solo_perf.carr,fsz);
ax2.YLabel.String       = 'State Accuracy [norm]';

% ACCURACY PER COHERENCE --> target-aligned
ax3                  	= axes('Position', [colmn(4) height(1) dim]); hold on
[dat,cond]              = transformData(solo_perf.acc_trg);
vl                      = violinplot(dat,cond);
vl                      = improveViolin(vl,ax3,col,solo_perf.carr,fsz);
ax3.YLabel.String       = 'Target Accuracy [norm]';

% ECCENTRICIY PER COHERENCE  --> state-aligned
ax4                  	= axes('Position', [colmn(3) height(2) dim]); hold on
[dat,cond]              = transformData(solo_perf.ecc_state);
vl                      = violinplot(dat,cond);
vl                      = improveViolin(vl,ax4,col,solo_perf.carr,fsz);
ax4.YLabel.String       = 'State Eccentricity [norm]';

% ECCENTRICIY PER COHERENCE --> target-aligned
ax5                  	= axes('Position', [colmn(4) height(2) dim]); hold on
[dat,cond]              = transformData(solo_perf.ecc_trg);
vl                      = violinplot(dat,cond);
vl                      = improveViolin(vl,ax5, col,solo_perf.carr,fsz);
ax5.YLabel.String       = 'Target Eccentricity [norm]';

% XCORR PER COHERENCE
ax6                  	= axes('Position', [colmn(1) height(3) dim]); hold on
ln                      = line([nLag nLag],[0 1],'Color',[0 0 0],'LineStyle', ':','LineWidth', lw);

for iCoh = 1:length(solo_perf.carr)
    cIdx                = solo_cr.coh == solo_perf.carr(iCoh);
    pline               = plot(mean(solo_cr.sxc(cIdx,:)),'Color', col(iCoh,:),'LineWidth',lw);
end

ax6.XLim             	= [1 300];
ax6.YLim              	= [0 round(max(mean(solo_cr.sxc(cIdx,:))),2)];
ax6.XTick             	= [1 150 300];
ax6.XTickLabel        	= {round([-150 0 150]*8.333)};
ax6.XLabel.String     	= 'Time relative to direction change [ms]';
ax6.YLabel.String    	= 'Avg. correlation coefficient';
ax6.FontSize          	= fsz;

% TARGET ECC-VS-ACC MATRIX
ax7                  	= axes('Position', [colmn(2) height(2) dim]); hold on

% Total juice: Target vs Cycle reward
ax8                  	= axes('Position', [colmn(1) height(2) dim]); hold on
reward_cumulative       = d.value(idx.reward);
% reward_cycle            = ???
% bar([reward_cumulative-reward_cycle reward_cycle])
% Moving average?

outcome               	= d.value(idx.outcome);
fixbreak_ratio        	= sum(strcmp(outcome, 'FixationBreak')) / length(cyc.cOn);

% Text displaying general information of training session
string_info = {['nCycles [#]: ' num2str(length(cyc.cOn))], ...
               ['nTargets [#]: ' num2str(length(solo_perf.trg_all_outc))], ...
               ['hir [%]: ' num2str(solo_perf.hir_pool)], ...
               ['juice [ml]: ' num2str(reward_cumulative{end})],...
               ['fixBreak [%]: ' num2str(fixbreak_ratio)]};

annotation('textbox', [0.35, 0.3, 0.2, 0.1], 'String', string_info, ...
           'FitBoxToText', 'on', 'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'FontSize',16);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Anahita - please add the functionality below %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GAZE HEAT MAP
ax9                  	= axes('Position', [colmn(3) height(3) dim]); hold on

% FIX BREAKS (absolute [#] and relative [%]) per coherence
ax10                  	= axes('Position', [colmn(4) height(3) dim]); hold on

% Save plots
print(f, [dest_dir_plot '/' fname(1:end-5)], '-r500', '-dpng');
print(f, [dest_dir_plot '/' fname(1:end-5)], '-r500', '-dsvg', '-painters');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Support Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = PHY_import_mwk2(fname, var_lst, write_file, cfg_pth)

if isfile([fname(1:end-5) '.h5']) && write_file == false                    % If file available...
    d       = MW_readH5([fname(1:end-5) '.h5']);                            % ...load .h5 file
else
    d     	= MW_readFile(fname, 'include', var_lst, '~typeOutcomeCheck','dotPositions'); % Import .mwk2 file
    
    if write_file
        disp('Save data structure...')
        MW_writeH5(d, [fname(1:end-5) '.h5'], 'replace', 'privateCFG', cfg_pth) % Save to .h5
        disp('Done!')
    end
end
end

function t = PHY_construct_table(subj, fname, d, cyc, idx, write_file)

% This function constructs a state-wise data table based on a previously
% imported .mwk2 data structure.
%
% Input:  	.subj           String, Subject ID
%           .fid            String, File identifier
%           .d              Structure, Contains raw data
%          	.cyc            Structure, Contains cycle information
%           .idx            Structure, Contains variable-wise indices
%           .write_file     Boolean, Indicates if table is to be saved
%
% Output:   .t              Table, State-wise data table containing
%                           stimulus parameters as well as behavioural
%                           responses
%
% Example:	t = CPR_construct_table('fxs','20210504',d, cyc, idx, true)
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
    'cyc_no', 'double'; ...                 % Trial/Block number
    'cyc_dur', 'double'; ...                % Trial/Block duration
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
    'fix_flag', 'double'; ...               % Fixation break flag
    'rdp_dot', 'cell'};                     % RDP dot position

% Initialise table
t       = table('Size',[5000, size(var,1)],...
    'VariableTypes',var(:,2),...
    'VariableNames',var(:,1));

%% Build table

fprintf('Building table...\n')

% Extract title information
fid          	= strsplit(fname,'_');

if size(fid,2) == 1
    setup       = 'physio4';
else
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
end

%%% Initial fixation check %%%
if ~strcmp(subj, 'agnt')
    fix_session = sum_duration_fix_breaks(d.value(idx.fixation), d.time(idx.fixation), d.time(end)-d.time(1)/1e3, d.time(end));
    if fix_session == true
        error('check_fixation')
    end
end

% Stimulus cycle (i.e. trial) loop
for iCyc = 1:length(cyc.cEnd)
    
    % Trial index
    cycIdx                  = [];
    cycIdx                  = d.time >= cyc.cOn(iCyc) & d.time <= cyc.cEnd(iCyc);
    
    % Extract dots
    tmp_dp                	= d.value(idx.RDP_dot);
    tmp_dp_ts           	= d.time(idx.RDP_dot);
    tmp_frme_ts             = d.time(idx.frame);
    
    % Stimulus trial data
    tmp_dir_ts              = getTrialData(d.time, cycIdx, idx.RDP_dir);                    % RDP_direction timestamps
    tmp_dir                 = getTrialData(d.value, cycIdx, idx.RDP_dir);                   % RDP_direction
    tmp_coh                 = getTrialData(d.value, cycIdx, idx.RDP_coh);                   % RDP coherence
    tmp_coh_ts            	= getTrialData(d.time, cycIdx, idx.RDP_coh);                    % RDP coherence timestamps
    
    % Remove first entries -> similar to last stimulus state at RDP_onset
    cyc.state_ts{iCyc}      = tmp_dir_ts(2:end);                                            % State ON duration timestamps
    cyc.state_dur{iCyc}    	= diff([tmp_dir_ts(2:end) cyc.cEnd(iCyc)]) ./ 1e3;              % Steady state duration for this trials [ms]
    cyc.dir{iCyc}           = tmp_dir(2:end);                                               % Trial RDP direction
    
    % Stimulus state loop
    for iSS = 1:length(cyc.state_dur{iCyc})
        
        % Steady state index
        if iSS < length(cyc.state_dur{iCyc})
            ssIdx           = d.time >= cyc.state_ts{iCyc}(iSS) & d.time < cyc.state_ts{iCyc}(iSS+1);
        else
            ssIdx           = d.time >= cyc.state_ts{iCyc}(iSS) & d.time <= cyc.cEnd(iCyc);
        end
        
        % Trial/State/Stimulus parameter
        cc                  = cc+1;                                                         % Steady state counter
        t.ID{cc}            = subj;                                                         % Subject ID
        t.date(cc)          = fid{1};                                                       % Experimental date
        t.exp{cc}           = fid{3};                                                       % Experimental condition
        t.setup{cc}         = setup;                                                        % Setup
        t.block(cc)         = str2double(fid{4}(end));                                      % Experimental block
        t.cyc_no(cc)        = iCyc;                                                         % Trial number
        t.cyc_dur(cc)       = (cyc.cEnd(iCyc) - cyc.cOn(iCyc)) ./ 1e3;                      % Trial duration [ms]
        t.ss_no(cc)         = cc;                                                           % Stimulus state counter
        t.ss_dur(cc)        = cyc.state_dur{iCyc}(iSS);                                     % Steady state duration
        t.rdp_coh(cc)       = tmp_coh(find(tmp_coh_ts <= cyc.state_ts{iCyc}(iSS),1,'last'));% Stimulus state coherence
        t.rdp_dir(cc)       = mod(cyc.dir{iCyc}(iSS),360);                                  % Stimulus direction
        
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
        
        % Sum duration of fixation breaks
        t.fix_flag(iSS)             = sum_duration_fix_breaks(d.value(ssIdx&idx.fixation), d.time(ssIdx&idx.fixation), t.ss_dur(cc), t.frme_ts{cc}(end));
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

% Remove fixation breaks - target outcomes only
outcome_ts(strcmp(outcome, 'FixationBreak'))    = [];
outcome(strcmp(outcome, 'FixationBreak'))       = [];

% Extract state and joystick parameter
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
    
%     tmp_outcome                     = outcome(find(outcome_ts >= trg_ts(iTarget),1,'first'));    % Extract target outcome
%     t.trg_hit{iState}               = [t.trg_hit{iState} , strcmp(tmp_outcome, 'hit')];         % Indicate if hit or miss
    t.trg_hit{iState}               = [t.trg_hit{iState} , strcmp(outcome(iTarget), 'hit')];         % Indicate if hit or miss
    
    smple                           = find(d.time(idx.JS_dir) <= trg_ts(iTarget), 1, 'last');  	% Extract last JS sample before target
    js_dev                          = rad2deg(circ_dist(deg2rad(js_dir_val{smple}),deg2rad(t.rdp_dir(iState)))); % Get circular distance to RDP direction
    js_acc                       	= abs(1 - abs(js_dev / 180));                               % Calculate accuracy
    
    t.trg_acc{iState}               = [t.trg_acc{iState} js_acc];                               % Add accuracy to vector
    t.trg_ecc{iState}               = [t.trg_ecc{iState} js_ecc_val{smple}];                    % Add eccentricity to vector
    
%     if strcmp(tmp_outcome, 'hit')
    if strcmp(outcome(iTarget), 'hit')
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

function out = sum_duration_fix_breaks(in_val, in_ts, duration, last_sample)
break_dur                   = [];

if ~isempty(in_val)    
    % Add duration of all fixation breaks
    for iBreak = find(cell2mat(in_val) == 0)
        if iBreak == length(in_val)
            break_stop = last_sample;
        else
            break_stop = in_ts(iBreak+1);
        end
        
        break_dur = [break_dur (break_stop - in_ts(iBreak))/1e3];               % Duration of fixation breaks
    end   
end

% Add fixation flag [boolean]
if ~isempty(in_val) && (sum(break_dur)/duration) > .1                           % If all breaks > 10% of state duration...
    out         = true;
else
    out         = false;
end
end

function [out, ps] = PHY_correlation_analysis(in, nLag, plotFlag)

% This function extracts stimulus and joystick data in a specified time
% window at the end of a steady state
%
% Input:  	.in             Table, Contains stimulus state information
%          	.nLag           Integer, Number of samples before direction
%                           changes (End of state)
%           .plotFlag       Boolean, Plot: yes - no
%
% Output:   .out            Structure, Contains averages and performance
%                           measures for each coherence level
%
% Example:	out = CPR_time_window(tbl,29,true)
%
% Known bugs:
%
% Feature wish list:
%
% Version history
%   1.0     (fxs 2021-05-04) Initial version.
%   1.1     (fxs 2021-05-25) Improved vectorisation procedure of data
%   1.2     (fxs 2021-11-25) Simplify code with frame-wise analysis

%% Check input

addpath /Users/fschneider/Documents/MATLAB/CircStat2012a

if nargin < 3 || isempty(plotFlag)
    plotFlag     	= true;
end

if nargin < 2 || isempty(nLag)
    nLag            = 150;
end

if nargin < 1
    error('Input missing')
end

%% Loop through blocks and correlate stimulus and behaviour

count                           = 0;

for iTrl = 1:size(in.frme_ts,2)
    
    % Sample-by-sample difference
    clear js_dff rdp_dff
    
    excl                        = isnan(in.rdp_coh{iTrl}) | isnan(in.rdp_dir{iTrl}) | isnan(in.js_dir{iTrl}) | isnan(in.frme_ts{iTrl});
    if ~isempty(excl)
        in.rdp_coh{iTrl}(excl)	= [];
        in.rdp_dir{iTrl}(excl) 	= [];
        in.js_dir{iTrl}(excl)  	= [];
        in.frme_ts{iTrl}(excl) 	= [];
    end
    
    rfr{iTrl}                   = median(in.refresh{iTrl});
    
    if isempty(in.frme_ts{iTrl})
        continue
    end
    
    for iSample = 1:size(in.frme_ts{iTrl},2)-1
        js_dff(iSample)     	= rad2deg(circ_dist(deg2rad(in.js_dir{iTrl}(iSample)),deg2rad(in.js_dir{iTrl}(iSample+1))));
        rdp_dff(iSample)     	= rad2deg(circ_dist(deg2rad(in.rdp_dir{iTrl}(iSample)),deg2rad(in.rdp_dir{iTrl}(iSample+1))));
    end
 
    js_dff                      = [0 js_dff];
    rdp_dff                     = [0 rdp_dff];
    
    % NANs coming from table-building or MWorks?!
    ex                          = isnan(js_dff) | isnan(rdp_dff);
    js_dff(ex)                  = [];
    rdp_dff(ex)                 = [];
    
    % Extract coherence chunks
    cindx                       = diff(in.rdp_coh{iTrl}) ~=0;
    cid                         = in.rdp_coh{iTrl}([true cindx]);
    cts                         = in.frme_ts{iTrl}([true cindx]);
            
    if length(cid) > 3 || sum(rdp_dff) == 0
        continue
    end
        
    for iCoh = 1:size(cts,2)
                
        % Build coherence index
        if iCoh < size(cts,2)
            cIdx                = in.frme_ts{iTrl} >= cts(iCoh) & in.frme_ts{iTrl} < cts(iCoh+1);
        else
            cIdx                = in.frme_ts{iTrl} >= cts(iCoh) & in.frme_ts{iTrl} <= in.frme_ts{iTrl}(end);
        end
        
        if sum(cIdx) < 10 %|| isempty(find(diff(rdp_dir{iTrl}(cIdx))))
            continue
        end
        
        % Correlation analysis
        count                   = count + 1;
        coh(count)              = cid(iCoh);                                                                % Coherence ID
        cc(count)               = circ_corrcc(deg2rad(in.js_dir{iTrl}(cIdx)),deg2rad(in.rdp_dir{iTrl}(cIdx)));	% Circular correlation between stimulus and joystick direction
        xc(count,:)             = xcorr(abs(js_dff(cIdx)),double(abs(rdp_dff(cIdx))),nLag,'normalized');   	% Cross-correlation between stimulus and joystick direction
        sxc(count,:)            = smoothdata(xc(count,:),'gaussian',20);                                    % Smooth correlation output with gaussian kernel
        maxR(count)           	= max(sxc(count,nLag:end));                                                 % Max cross-correlation coefficient
        posPk(count)        	= find(sxc(count,nLag:end) == maxR(count));                                 % Peak position of cross-correlation
        
        fs(count)               = rfr{iTrl};
        lg(count)             	= (find(sxc(count,nLag:end) == max(sxc(count,nLag:end))) * fs(count)) / 1e3; % Calculate lag using sample rate

        try
            auPk(count,:)     	= trapz(sxc(count,posPk(count)-10:posPk(count)+10));                        % Area under cross-correlation peak
        catch
            auPk(count,:)     	= nan;
        end
        
        if plotFlag
            % Plot trial-wise cross-correlation
            ps              	= plot(sxc(count,nLag:end),'Color',[.5 .5 .5 .1], 'Marker','none', 'LineWidth',2);
        else 
            ps                  = [];
        end
    end
end

out.coh                         = coh;
out.cc                          = cc;
out.xc                          = xc;
out.sxc                         = sxc;
out.maxR                        = maxR;
out.posPk                       = posPk;
out.auPk                        = auPk;
out.lag                       	= lg;

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
    tmp_trg_ts_state{j} 	= (tmp_trg_ts{j} - tmp_frame1_ts(j)) ./ 1e3;  
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
    trgIdx                  = (t1_ts-f1_ts) >= 1e6;
    rdp_dir                 = in.rdp_dir(cIdx & in.trg_shown & trgIdx);
    js_dir                  = in.js_dir(cIdx & in.trg_shown & trgIdx);
    js_ecc_tmp              = in.js_ecc(cIdx & in.trg_shown & trgIdx);
    frmes                   = in.frme_ts(cIdx & in.trg_shown & trgIdx);
    trg1_ts                 = t1_ts(cIdx & in.trg_shown & trgIdx);
    
    clear js_acc js_ecc
    for iTrg = 1:length(rdp_dir)
        clear js_dev
        smpl_idx            = find(frmes{iTrg} < trg1_ts(iTrg),1,'last')-nSample : find(frmes{iTrg} < trg1_ts(iTrg),1,'last');
        js_dev              = rad2deg(circ_dist(deg2rad(js_dir{iTrg}(smpl_idx)),deg2rad(rdp_dir(iTrg))));  % Minimum RDP-Joystick difference
        js_acc(iTrg)      = nanmean(abs(1 - abs(js_dev) / 180));         	% Joystick accuracy
        js_ecc(iTrg)      = nanmean(js_ecc_tmp{iTrg}(smpl_idx));
    end
    
    out.macc_trg(iCoh)      = nanmedian(js_acc);
    out.mecc_trg(iCoh)      = nanmedian(js_ecc);
    out.acc_trg{iCoh}     	= js_acc;
    out.ecc_trg{iCoh}       = js_ecc;
    out.outc_trg{iCoh}      = cell2mat(t1_outc(cIdx & in.trg_shown & trgIdx));
    out.carr(iCoh)         	= snr(iCoh);
end
end

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

function vl = improveViolin(vl,ax,col_map,snr,fsz)

for iV=1:length(vl)
    if isempty(vl(iV).ScatterPlot)
        continue
    else
        vl(iV).BoxWidth                     = .03;
        vl(iV).ViolinColor{1}               = col_map(iV,:);
        vl(iV).ViolinAlpha{1}               = .5;
        vl(iV).ScatterPlot.MarkerFaceColor  = [.25 .25 .25];
        vl(iV).ScatterPlot.MarkerEdgeColor  = 'none';
        vl(iV).ScatterPlot.SizeData         = 25;
    end
end

ax.XTickLabel    	= {snr};
ax.XLabel.String   	= 'Coherence';
ax.YLim           	= [0 1];
ax.FontSize         = fsz;
box off
end

function [dat, cond]= transformData(in)
    dat = []; cond = [];
    
    for iCoh = 1:length(in)
        if size(in{iCoh},1) == 1
            dat            	= [dat; in{iCoh}'];
        else
            dat          	= [dat; in{iCoh}];
        end
        cond             	= [cond; repmat(iCoh,length(in{iCoh}),1)];
    end
end

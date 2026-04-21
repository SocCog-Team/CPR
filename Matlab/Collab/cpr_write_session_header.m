function cpr_write_session_header(phy, output_dir)
% CPR_WRITE_SESSION_HEADER  Write README.txt and per-cycle JSON header files.
%
%   CPR_WRITE_SESSION_HEADER(PHY, OUTPUT_DIR) generates two types of files
%   that document a CPR recording session:
%
%       README.txt
%           Human-readable overview of the full session: recording metadata,
%           task description, folder structure, and schema documentation for
%           all output file types. Intended as the entry point for anyone
%           exploring the dataset.
%
%       cycle_headers/<rec>_<block>_cycle<NNN>_header.json
%           One JSON file per experimental cycle containing all metadata
%           relevant to that cycle: condition, subjects, stimulus parameters,
%           timing, and a behavioural summary (hit/miss/fixation-break counts).
%
%   This function is designed to be called once per recording session, before
%   or after cpr_write_cycle_frames and cpr_write_unit_spikes.
%
% -------------------------------------------------------------------------
%   INPUTS
%       phy         Master experiment struct. Required fields:
%
%                   phy.exp.date            Recording date string 'YYYYMMDD'
%                   phy.exp.monkey          Monkey subject ID
%                   phy.exp.experimenter    Experimenter / human subject ID
%                   phy.exp.block           Block label, e.g. 'block1'
%                   phy.exp.setup           Setup identifier, e.g. 'phy4'
%                   phy.exp.rec_num         Recording ID, e.g. 'rec068'
%
%                   phy.exp.human           (optional) Human subject ID.
%                                           Falls back to phy.exp.experimenter.
%
%                   phy.stim.cpr_cyle       [nCycles × 2] cycle start/end (s)
%                   phy.stim.cpr_solo       {1 × nCycles} logical cells
%                   phy.stim.cpr_dyad       {1 × nCycles} logical cells
%                   phy.stim.task           {1 × nCycles} task label strings
%                   phy.stim.rdp_center_xy  {1 × nCycles} [1×2] position (deg)
%                   phy.stim.rdp_dir        {1 × nCycles} direction vectors
%                   phy.stim.rdp_coh        {1 × nCycles} coherence vectors
%                   phy.stim.frme_ts        {1 × nCycles} frame timestamps
%                   phy.stim.feedback_ts    {1 × nCycles} feedback timestamps
%                   phy.stim.outcome        {1 × nCycles} outcome string cells
%                   phy.stim.reward_cum     {1 × nCycles} cumulative reward
%
%       output_dir  Root output folder path. Created if it does not exist.
%
% -------------------------------------------------------------------------
%   OUTPUT STRUCTURE
%       <output_dir>/
%           README.txt
%           cycle_headers/
%               <rec>_<block>_cycle001_header.json
%               <rec>_<block>_cycle002_header.json
%               ...
%
% -------------------------------------------------------------------------
%   REQUIREMENTS
%       MATLAB R2021a or later (jsonencode 'PrettyPrint' option).
%       For earlier versions, replace jsonencode(..., 'PrettyPrint', true)
%       with jsonencode(...) — output will be valid but not formatted.
%
% -------------------------------------------------------------------------
%   USAGE
%       cpr_write_session_header(phy, fullfile('export', phy.exp.rec_num))
%
%   See also: cpr_write_cycle_frames, cpr_write_unit_spikes
%
% -------------------------------------------------------------------------
%   Created    : 2026-04-16
%   AI Source  : Generated with Claude (Anthropic) — https://claude.ai
%   Reference  : Schneider et al. (2025) eLife https://doi.org/10.7554/eLife.101021.2

% =========================================================================
% 0.  SETUP
% =========================================================================

header_dir = fullfile(output_dir, 'cycle_headers');
ensure_dir(output_dir);
ensure_dir(header_dir);

% Use cpr_cyle row count as the authoritative cycle count.
% cpr_cyle is [nCycles x 2] — one row per cycle with [t_start, t_end].
n_cycles = size(phy.stim.cpr_cyle, 1);

% =========================================================================
% 1.  README
% =========================================================================

write_readme(phy, output_dir, n_cycles);

% =========================================================================
% 2.  PER-CYCLE JSON HEADERS
% =========================================================================

for iCycle = 1:n_cycles
    write_cycle_header_json(phy, iCycle, output_dir, header_dir);
end

fprintf('[cpr_write_session_header] README + %d cycle header(s) written to: %s\n', ...
    n_cycles, output_dir);

end % main function


% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function write_readme(phy, output_dir, n_cycles)
% WRITE_README  Write the human-readable README.txt for the dataset root.
%   Documents the session, folder layout, and all file schemas so that a
%   new collaborator can understand the dataset without needing MATLAB.

e        = phy.exp;
rec      = e.rec_num;
blk      = e.block;
gen_date = datestr(now, 'yyyy-mm-dd');

readme_path = fullfile(output_dir, 'README.txt');
fid = open_file(readme_path);

ln  = @(s) fprintf(fid, '%s\n', s);    % print line helper
sep = @()  ln(repmat('=', 1, 62));      % section separator

sep();
ln(' CPR NEURAL AND BEHAVIOURAL DATASET');
sep();
fprintf(fid, ' Recording   : %s\n', rec);
fprintf(fid, ' Date        : %s\n', e.date);
fprintf(fid, ' Monkey      : %s\n', e.monkey);
fprintf(fid, ' Experimenter: %s\n', e.experimenter);
fprintf(fid, ' Setup       : %s\n', e.setup);
fprintf(fid, ' Block       : %s\n', blk);
fprintf(fid, ' N cycles    : %d\n', n_cycles);
fprintf(fid, ' Generated   : %s\n', gen_date);
fprintf(fid, ' AI source   : Claude (Anthropic) — https://claude.ai\n');
ln('');

% --- Overview ------------------------------------------------------------
ln('OVERVIEW');
ln('--------');
ln(' This dataset contains neural spike times and behavioural responses from');
ln(' a Continuous Perceptual Report (CPR) task. Subjects report the perceived');
ln(' direction of a Random Dot Pattern (RDP) in real time using a joystick.');
ln(' Joystick angle = perceived direction; joystick tilt = confidence.');
ln(' Sessions run in either solo (monkey only) or dyadic (monkey + human)');
ln(' conditions. Neural data are extracellular single-unit spike times.');
ln('');
ln(' Reference: Schneider et al. (2025) eLife');
ln('            https://doi.org/10.7554/eLife.101021.2');
ln('');

% --- Folder structure ----------------------------------------------------
ln('FOLDER STRUCTURE');
ln('----------------');
ln('  README.txt                                  This file');
ln('  cycle_headers/');
fprintf(fid, '      %s_%s_cycle001_header.json     Metadata per cycle\n', rec, blk);
ln('      ...');
ln('  cycle_frames/');
fprintf(fid, '      %s_%s_cycle001_frames.csv      Frame-aligned data per cycle\n', rec, blk);
ln('      ...');
ln('  unit_spikes/');
fprintf(fid, '      %s_ch001_neg_unit1_range1_spikes.csv   Spike times per unit range\n', rec);
ln('      ...');
ln('');

% --- Cycle header schema -------------------------------------------------
ln('CYCLE HEADER SCHEMA  (cycle_headers/*.json)');
ln('--------------------------------------------');
ln('  recording           Recording identifier');
ln('  date                Recording date (YYYYMMDD)');
ln('  experimenter        Experimenter ID');
ln('  setup               Setup identifier');
ln('  block               Block label');
ln('  cycle               Cycle number (1-based)');
ln('  task                Task type string');
ln('  condition           "solo" or "dyad"');
ln('  cpr_solo            Boolean: true if solo condition');
ln('  cpr_dyad            Boolean: true if dyadic condition');
ln('  subject_monkey      Monkey subject ID');
ln('  subject_human       Human subject ID');
ln('  rdp_center_x_deg    RDP centre x-position (visual degrees)');
ln('  rdp_center_y_deg    RDP centre y-position (visual degrees)');
ln('  cycle_start_s       Absolute cycle onset timestamp (seconds)');
ln('  cycle_end_s         Absolute cycle end timestamp (seconds)');
ln('  cycle_duration_s    Cycle duration (seconds)');
ln('  n_frames            Number of display frames in cycle');
ln('  n_direction_changes Number of RDP direction change events');
ln('  n_coherence_changes Number of RDP coherence change events');
ln('  n_feedback_events   Total feedback events');
ln('  n_hits              Number of hit outcomes');
ln('  n_misses            Number of miss outcomes');
ln('  n_fixation_breaks   Number of FixationBreak events');
ln('  fixation_break      Boolean: true if cycle ended with fixation break');
ln('  reward_total        Cumulative reward at end of cycle');
ln('  frame_rate_hz       Display refresh rate (Hz)');
ln('  frames_file         Relative path to the cycle frame CSV');
ln('');

% --- Frame CSV schema ----------------------------------------------------
ln('FRAME CSV SCHEMA  (cycle_frames/*.csv)');
ln('---------------------------------------');
ln('  Each row is one display frame at 120 Hz.');
ln('');
ln('  frame_idx              0-based integer frame counter');
ln('  time_s                 Time from cycle onset (seconds)');
ln('  stim_direction_deg     RDP direction [0, 360) degrees');
ln('                           Step-held between change events; not interpolated');
ln('  stim_coherence         RDP coherence fraction [0.0, 1.0]');
ln('                           Step-held between change events; not interpolated');
ln('  resp_dir_monkey_deg    Monkey joystick direction [0, 360) degrees');
ln('                           Circular linear interpolation to frame grid');
ln('  resp_conf_monkey       Monkey joystick tilt magnitude [0.0, 1.0]');
ln('                           Scalar linear interpolation to frame grid');
ln('  resp_dir_human_deg     Human joystick direction [dyad only; NaN in solo]');
ln('  resp_conf_human        Human joystick tilt      [dyad only; NaN in solo]');
ln('  feedback_outcome_code  Priority-encoded outcome at feedback frames:');
ln('                            1  = hit');
ln('                            0  = miss');
ln('                           -1  = FixationBreak');
ln('                           NaN = no feedback event on this frame');
ln('  feedback_outcome_str   All outcomes at this frame, pipe-separated.');
ln('                           e.g. "hit", "miss", "hit|FixationBreak"');
ln('                           Empty string on non-feedback frames.');
ln('  reward_ind             Reward earned at this specific feedback event.');
ln('                           miss -> 0; hit -> positive amount.');
ln('                           NaN on non-feedback frames and FixationBreak-only events.');
ln('  reward_cycle_cum       Running reward total within the current cycle.');
ln('                           Resets to 0 at the start of each cycle.');
ln('                           NaN on non-feedback frames.');
ln('  reward_session_cum     Session-wide cumulative reward at this feedback event.');
ln('                           Raw values from phy.stim.reward_cum, never zeroed.');
ln('                           NaN on non-feedback frames.');
ln('');
ln('  OUTCOME PRIORITY (when multiple outcomes share a frame):');
ln('    hit (1) > miss (0) > FixationBreak (-1)');
ln('    All outcomes are listed in feedback_outcome_str regardless of priority.');
ln('');
ln('  NOTE ON FIXATION BREAKS:');
ln('    A FixationBreak terminates the cycle prematurely. It may appear without');
ln('    a matching reward timestamp (it is a cycle-end signal, not a reward event).');
ln('    In that case it is snapped to the last frame with NaN for both reward columns.');
ln('');

% --- Spike CSV schema ----------------------------------------------------
ln('SPIKE CSV SCHEMA  (unit_spikes/*.csv)');
ln('--------------------------------------');
ln('  Long format: one row per spike.');
ln('  File name encodes channel, unit, and recording range.');
ln('');
ln('  cycle_id        Cycle number (1-based, matches frame CSV index)');
ln('  spike_time_us   Spike time in microseconds, relative to cycle onset');
ln('');
ln('  NOTES:');
ln('  - Only cycles in which the unit was present (cyc_id mask) are included.');
ln('  - Cycles with zero spikes produce no rows.');
ln('  - spike_time_us / 1e6 = spike_time_s, comparable to time_s in frame CSV.');
ln('  - Only significantly stimulus-driven units are exported');
ln('    (signrank test; see inclusion_flag in phy.brain.CPR.spks.include).');
ln('  - "range" in the filename denotes the span of cycles over which the');
ln('    unit was tracked as a stable cluster (e.g. range1 = cycles 1-100,');
ln('    range2 = cycles 132-177). Same channel+unit across ranges is the');
ln('    same sorted cluster recorded at different times within the session.');
ln('');

sep();
fclose(fid);

end % write_readme


% -------------------------------------------------------------------------

function write_cycle_header_json(phy, iCycle, output_dir, header_dir)
% WRITE_CYCLE_HEADER_JSON  Assemble and write a single cycle header JSON.
%   Collects all metadata for cycle ICYCLE: timing, stimulus summary,
%   subject IDs, and a behavioural outcome count.

e   = phy.exp;
stm = phy.stim;

% --- Condition -----------------------------------------------------------
is_solo = stm.cpr_solo{iCycle};
is_dyad = stm.cpr_dyad{iCycle};

if     is_solo && ~is_dyad,  condition = 'solo';
elseif is_dyad && ~is_solo,  condition = 'dyad';
else,                         condition = 'unknown';
    warning('cpr_write_session_header:ambiguousCondition', ...
        'Cycle %d: cpr_solo and cpr_dyad are both %d. Condition set to ''unknown''.', ...
        iCycle, is_solo);
end

% --- Subject IDs ---------------------------------------------------------
% Subject_human is always recorded; condition clarifies active participation.
monkey_id = e.monkey;
human_id  = get_field_or_default(e, 'human', e.experimenter);

% --- Task type -----------------------------------------------------------
if iscell(stm.task) && iCycle <= numel(stm.task)
    task_str = char(stm.task{iCycle});
else
    task_str = 'CPR';
end

% --- RDP screen position -------------------------------------------------
rdp_xy = double(stm.rdp_center_xy{iCycle});   % [1 x 2] in visual degrees

% --- Cycle timing --------------------------------------------------------
cycle_ts      = double(stm.cpr_cyle(iCycle,:));
cycle_start_s = cycle_ts(1);
cycle_end_s   = cycle_ts(2);
cycle_dur_s   = cycle_end_s - cycle_start_s;

% --- Frame count ---------------------------------------------------------
n_frames = numel(stm.frme_ts{iCycle});

% --- Stimulus event counts -----------------------------------------------
n_dir_changes = numel(stm.rdp_dir{iCycle});
n_coh_changes = numel(stm.rdp_coh{iCycle});

% --- Outcome summary -----------------------------------------------------
outcome_cell = stm.outcome{iCycle};   % {1 x nEvents} cell of strings
n_hits       = sum(strcmpi(outcome_cell, 'hit'));
n_misses     = sum(strcmpi(outcome_cell, 'miss'));
n_fixbreaks  = sum(strcmpi(outcome_cell, 'FixationBreak'));
n_feedback   = numel(outcome_cell);
fixbreak_flag = logical(n_fixbreaks > 0);

% --- Cycle-relative reward -----------------------------------------------
% reward_cum{iCycle} contains session-wide cumulative values. Subtracting
% the first value of the cycle's own reward_cum gives the reward earned
% within the cycle. The first value is the baseline because it is recorded
% at (or before) the first feedback event, before any reward is earned.
% An empty reward_cum means no feedback occurred this cycle: reward = 0.
reward_cv_session = double(stm.reward_cum{iCycle});
if ~isempty(reward_cv_session) && ~all(isnan(reward_cv_session))
    reward_total = reward_cv_session(end) - reward_cv_session(1);
else
    reward_total = 0;
end

% --- Relative path to the companion frame CSV ----------------------------
frames_file = fullfile('cycle_frames', ...
    sprintf('%s_%s_cycle%03d_frames.csv', e.rec_num, e.block, iCycle));

% --- Assemble header struct (field order = JSON key order) ---------------
h = struct();
h.recording           = e.rec_num;
h.date                = e.date;
h.experimenter        = e.experimenter;
h.setup               = e.setup;
h.block               = e.block;
h.cycle               = iCycle;
h.task                = task_str;
h.condition           = condition;
h.cpr_solo            = logical(is_solo);
h.cpr_dyad            = logical(is_dyad);
h.subject_monkey      = monkey_id;
h.subject_human       = human_id;
h.rdp_center_x_deg    = rdp_xy(1);
h.rdp_center_y_deg    = rdp_xy(2);
h.cycle_start_s       = cycle_start_s;
h.cycle_end_s         = cycle_end_s;
h.cycle_duration_s    = cycle_dur_s;
h.n_frames            = n_frames;
h.n_direction_changes = n_dir_changes;
h.n_coherence_changes = n_coh_changes;
h.n_feedback_events   = n_feedback;
h.n_hits              = n_hits;
h.n_misses            = n_misses;
h.n_fixation_breaks   = n_fixbreaks;
h.fixation_break      = fixbreak_flag;
h.reward_total        = reward_total;
h.frame_rate_hz       = 120;
h.frames_file         = frames_file;

% --- Write JSON ----------------------------------------------------------
json_name = sprintf('%s_%s_cycle%03d_header.json', e.rec_num, e.block, iCycle);
json_path = fullfile(header_dir, json_name);
fid = open_file(json_path);
fprintf(fid, '%s\n', jsonencode(h, 'PrettyPrint', true));
fclose(fid);

end % write_cycle_header_json


% -------------------------------------------------------------------------

function val = get_field_or_default(s, field, default)
% GET_FIELD_OR_DEFAULT  Return s.field if it exists, otherwise default.
if isfield(s, field)
    val = s.(field);
else
    val = default;
end
end % get_field_or_default


% -------------------------------------------------------------------------

function fid = open_file(filepath)
% OPEN_FILE  Open a file for writing, raising a clear error on failure.
fid = fopen(filepath, 'w');
if fid == -1
    error('cpr_write_session_header:fileOpenFailed', ...
        'Could not open file for writing: %s', filepath);
end
end % open_file


% -------------------------------------------------------------------------

function ensure_dir(d)
% ENSURE_DIR  Create directory D if it does not already exist.
if ~exist(d, 'dir')
    mkdir(d);
end
end % ensure_dir


% -------------------------------------------------------------------------


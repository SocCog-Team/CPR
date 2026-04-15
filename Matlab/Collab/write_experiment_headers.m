function write_experiment_headers(phy, output_dir)
% WRITE_EXPERIMENT_HEADERS  Write README and per-cycle JSON headers for CPR data sharing.
%
%   WRITE_EXPERIMENT_HEADERS(PHY, OUTPUT_DIR) generates two file types:
%
%       README.txt              - General experiment & task description.
%                                 Written once per recording.
%
%       header_cycle{N}.json    - Cycle-specific metadata including condition,
%                                 subject IDs and the associated CSV filename.
%                                 Written for each cycle in phy.stim.cpr_solo.
%
% -------------------------------------------------------------------------
%   INPUTS
%       phy         Experiment struct. Required fields:
%
%                   phy.exp.date            Recording date, string 'YYYYMMDD'
%                   phy.exp.monkey          Subject ID
%                   phy.exp.block           Block identifier, e.g. 'block1'
%                   phy.exp.setup           Recording setup, e.g. 'phy4'
%                   phy.exp.rec_num         Recording number, e.g. 'rec068'
%                   phy.exp.experimenter    Experimenter ID, e.g. 'fxs'
%
%                   phy.stim.cpr_solo       {nCycles x 1} logical cell array
%                   phy.stim.cpr_dyad       {nCycles x 1} logical cell array
%
%                   phy.exp.subject_ids     (optional) struct with fields:
%                                               .p1  ID of partner 1
%                                               .p2  ID of partner 2
%                                           If absent, phy.exp.monkey is used
%                                           for p1 and p2 is set to null.
%
%       output_dir  Path to destination folder. Created if it does not exist.
%
% -------------------------------------------------------------------------
%   REQUIREMENTS
%       jsonencode with 'PrettyPrint' requires MATLAB R2021a or later.
%
% -------------------------------------------------------------------------
%   EXAMPLE
%       write_experiment_headers(phy, fullfile('data', phy.exp.rec_num))
%
% -------------------------------------------------------------------------
%   OUTPUT FILES
%       README.txt
%       header_cycle01.json, header_cycle02.json, ...
%
%   See also: jsonencode, jsondecode
%
%   Reference: https://doi.org/10.7554/eLife.101021.2

% --- output directory -------------------------------------------------
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

nCycles = numel(phy.stim.cpr_solo);

% --- write general README (once per recording) ------------------------
write_readme(phy, output_dir, nCycles);

% --- write per-cycle header (once per cycle) --------------------------
for iCycle = 1:nCycles
    write_cycle_header(phy, iCycle, output_dir);
end

fprintf('[write_experiment_headers] Done. %d cycle header(s) written to: %s\n', ...
    nCycles, output_dir);

end % main function


% =========================================================================
function write_readme(phy, output_dir, nCycles)
% Writes README.txt describing the full recording session and CSV schema.

readme_path = fullfile(output_dir, 'README.txt');
fid = open_file_for_writing(readme_path);

e = phy.exp;

fprintf(fid, '=======================================================\n');
fprintf(fid, ' CONTINUOUS PERCEPTUAL REPORT (CPR) — EXPERIMENT README\n');
fprintf(fid, '=======================================================\n\n');

fprintf(fid, 'EXPERIMENT\n');
fprintf(fid, '----------\n');
fprintf(fid, '  Date          : %s\n', e.date);
fprintf(fid, '  Subject       : %s\n', e.monkey);
fprintf(fid, '  Block         : %s\n', e.block);
fprintf(fid, '  Setup         : %s\n', e.setup);
fprintf(fid, '  Recording     : %s\n', e.rec_num);
fprintf(fid, '  Experimenter  : %s\n', e.experimenter);
fprintf(fid, '  N cycles      : %d\n\n', nCycles);

fprintf(fid, 'TASK\n');
fprintf(fid, '----\n');
fprintf(fid, '  Name          : Continuous Perceptual Report (CPR)\n');
fprintf(fid, '  Stimulus      : Random Dot Pattern (RDP)\n');
fprintf(fid, '  Response      : Joystick\n');
fprintf(fid, '                    direction = perceived motion angle (degrees)\n');
fprintf(fid, '                    tilt      = confidence (normalised [0, 1])\n');
fprintf(fid, '  Frame rate    : 120 Hz\n');
fprintf(fid, '  Conditions    : solo (cpr_solo) / dyadic (cpr_dyad)\n');
fprintf(fid, '  Reference     : https://doi.org/10.7554/eLife.101021.2\n\n');

fprintf(fid, 'FILE STRUCTURE\n');
fprintf(fid, '--------------\n');
fprintf(fid, '  README.txt                   : this file\n');
fprintf(fid, '  header_cycle{NN}.json        : cycle-specific metadata\n');
fprintf(fid, '  %s_%s_cycle{NN}.csv          : frame-wise time-series data\n\n', ...
    e.rec_num, e.block);

fprintf(fid, 'CSV COLUMN SCHEMA\n');
fprintf(fid, '-----------------\n');
fprintf(fid, '  frame_idx           Integer frame index (0-based)\n');
fprintf(fid, '  time_s              Time from cycle onset in seconds\n');
fprintf(fid, '  stim_direction_deg  RDP motion direction in degrees [0, 360)\n');
fprintf(fid, '  stim_coherence      RDP coherence fraction [0.0, 1.0]\n');
fprintf(fid, '  resp_dir_monk_deg   Monkey joystick direction in degrees [0, 360)\n');
fprintf(fid, '  resp_conf_monk      Monkey joystick tilt magnitude [0.0, 1.0]\n');
fprintf(fid, '  resp_dir_hum_deg    Human  joystick direction [dyad only; NaN in solo]\n');
fprintf(fid, '  resp_conf_hum       Human  joystick tilt magnitude [dyad only; NaN in solo]\n');
fprintf(fid, '  feedback_outcome    Reward feedback event: 1=hit, 0=miss, NaN=no event\n');
fprintf(fid, '  reward_cum          Cumulative reward at feedback event (NaN on other frames)\n\n');

fprintf(fid, '=======================================================\n');
fclose(fid);

end % write_readme


% =========================================================================
function write_cycle_header(phy, iCycle, output_dir)
% Writes header_cycle{NN}.json for a single cycle.

e        = phy.exp;
is_solo  = phy.stim.cpr_solo{iCycle};
is_dyad  = phy.stim.cpr_dyad{iCycle};

% --- condition label --------------------------------------------------
if     is_solo && ~is_dyad,  condition = 'solo';
elseif is_dyad && ~is_solo,  condition = 'dyad';
else
    condition = 'unknown';
    warning('write_experiment_headers:ambiguousCondition', ...
        'Cycle %d: cpr_solo=%d and cpr_dyad=%d — condition set to ''unknown''.', ...
        iCycle, is_solo, is_dyad);
end

% --- subject IDs ------------------------------------------------------
if isfield(e, 'subject_ids') && isstruct(e.subject_ids)
    p1_id = e.subject_ids.p1;
    p2_id = e.subject_ids.p2;
else
    p1_id = e.monkey;
    p2_id = [];          % encodes as null in JSON
end

% In solo condition p2 is always null regardless of subject_ids
if is_solo && ~is_dyad
    p2_id = [];
end

% --- RDP screen position ----------------------------------------------
% rdp_center_xy: [x, y] position of the RDP centre in visual degrees
rdp_xy = double(phy.stim.rdp_center_xy{iCycle});   % [1 x 2] single -> double

% --- assemble header struct -------------------------------------------
header.recording        = e.rec_num;
header.date             = e.date;
header.experimenter     = e.experimenter;
header.setup            = e.setup;
header.block            = e.block;
header.cycle            = iCycle;
header.condition        = condition;
header.cpr_solo         = logical(is_solo);
header.cpr_dyad         = logical(is_dyad);
header.subject_p1       = p1_id;
header.subject_p2       = p2_id;
header.rdp_center_x_deg = rdp_xy(1);   % horizontal position in visual degrees
header.rdp_center_y_deg = rdp_xy(2);   % vertical position in visual degrees
header.frame_rate_hz    = 120;
header.data_file        = sprintf('%s_%s_cycle%02d.csv', e.rec_num, e.block, iCycle);

% --- write JSON -------------------------------------------------------
header_path = fullfile(output_dir, sprintf('header_cycle%02d.json', iCycle));
fid = open_file_for_writing(header_path);
fprintf(fid, '%s\n', jsonencode(header, 'PrettyPrint', true));
fclose(fid);

end % write_cycle_header


% =========================================================================
function fid = open_file_for_writing(filepath)
% Helper: open file, error on failure.

fid = fopen(filepath, 'w');
if fid == -1
    error('write_experiment_headers:fileOpenFailed', ...
        'Could not open file for writing: %s', filepath);
end

end % open_file_for_writing

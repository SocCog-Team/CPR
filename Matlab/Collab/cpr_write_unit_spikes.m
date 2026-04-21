function cpr_write_unit_spikes(phy, output_dir)
% CPR_WRITE_UNIT_SPIKES  Export spike times for significantly driven units.
%
%   CPR_WRITE_UNIT_SPIKES(PHY, OUTPUT_DIR) reads the spike inclusion table
%   and spike time structs, and writes one long-format CSV per unit range
%   to:
%       <output_dir>/unit_spikes/<rec>_<channel>_<unit>_<range>_spikes.csv
%
%   Each CSV contains all spike times for the cycles in which that unit
%   range was present, as specified by the cyc_id mask in the include table.
%
% -------------------------------------------------------------------------
%   INCLUSION CRITERIA
%
%       Only rows where inclusion_flag == true are exported.
%       Inclusion is based on a Wilcoxon signed-rank test comparing firing
%       rate during vs. before stimulus onset (see phy.brain.CPR.spks.include).
%
% -------------------------------------------------------------------------
%   UNIT RANGES
%
%       A physical neuron may be tracked as a stable spike cluster across
%       different, non-overlapping spans of cycles. Each span is stored as
%       a separate row in the include table with a range suffix in the
%       unit_ID (e.g. unit1_range1, unit1_range2). Each range receives its
%       own output file. This preserves the information that the unit was
%       re-identified or re-sorted across session blocks.
%
%       unit_ID format: <channel>_<unit_field>_<range>
%           e.g.  "ch001_neg_unit1_range1"
%           -> struct path: phy.brain.CPR.spks.ch001_neg.unit1{iCycle}
%           -> output file: <rec>_ch001_neg_unit1_range1_spikes.csv
%
% -------------------------------------------------------------------------
%   OUTPUT FORMAT  (long-format CSV with comment header)
%
%       Columns:
%           cycle_id        Cycle number (1-based; matches frame CSV index)
%           spike_time_us   Spike time in microseconds, relative to cycle onset
%
%       Each row is one spike. Long format was chosen over per-cycle files
%       for compactness — all data for a unit range is in one place.
%
%       Cycles listed in cyc_id but containing zero spikes produce no rows.
%       The absence of rows for a cycle indicates silence, not missing data.
%
%       The file begins with comment lines (prefixed '#') documenting unit
%       metadata. These are ignored by standard parsers:
%           Python:  pd.read_csv(f, comment='#')
%           MATLAB:  readtable(f, 'CommentStyle', '#')
%           R:       read.csv(f, comment.char='#')
%
% -------------------------------------------------------------------------
%   SPIKE TIME UNITS
%
%       spike_time_us is in microseconds relative to the cycle onset time
%       (phy.stim.cpr_cyle(iCycle, 1)). To convert to seconds for alignment
%       with the time_s column in the frame CSV:
%           spike_time_s = spike_time_us / 1e6
%
% -------------------------------------------------------------------------
%   INPUTS
%       phy         Master experiment struct. Required fields:
%                   phy.exp.rec_num, phy.exp.block
%                   phy.brain.CPR.spks.include  — 301×8 inclusion table
%                   phy.brain.CPR.spks.<ch>.<unit> — {1×nCycles} spike cells
%
%       output_dir  Root output folder. Created if absent.
%
% -------------------------------------------------------------------------
%   REQUIREMENTS
%       MATLAB R2016b or later (string arrays).
%
% -------------------------------------------------------------------------
%   USAGE
%       cpr_write_unit_spikes(phy, fullfile('export', phy.exp.rec_num))
%
%   See also: cpr_write_session_header, cpr_write_cycle_frames
%
% -------------------------------------------------------------------------
%   Created    : 2026-04-16
%   AI Source  : Generated with Claude (Anthropic) — https://claude.ai
%   Reference  : Schneider et al. (2025) eLife https://doi.org/10.7554/eLife.101021.2

% =========================================================================
% 0.  SETUP
% =========================================================================

spikes_dir = fullfile(output_dir, 'unit_spikes');
if ~exist(output_dir,  'dir'), mkdir(output_dir);  end
if ~exist(spikes_dir,  'dir'), mkdir(spikes_dir);  end

tbl  = phy.brain.CPR.spks.include;     % inclusion table (nUnits x 8)
spks = phy.brain.CPR.spks;             % struct containing channel substruct
rec  = phy.exp.rec_num;

% Identify rows that pass the stimulus-driven inclusion criterion
included_rows = find(tbl.inclusion_flag);

fprintf('[cpr_write_unit_spikes] %d / %d unit ranges pass inclusion.\n', ...
    numel(included_rows), height(tbl));

% =========================================================================
% 1.  ITERATE OVER INCLUDED UNIT RANGES
%
%   Each row in the include table represents one (unit, range) pair.
%   We process them individually — no range merging — so each gets its
%   own output file.
% =========================================================================

n_written = 0;

for iRow = included_rows(:)'

    % --- Parse unit_ID ---------------------------------------------------
    % Format: "<channel>_<unit_field>_<range>"
    % e.g.   "ch001_neg_unit1_range1"
    %          channel    = 'ch001_neg'
    %          unit_field = 'unit1'
    %          range_str  = 'range1'
    unit_id = char(tbl.unit_ID(iRow));

    [channel, unit_field, range_str, parse_ok] = parse_unit_id(unit_id);

    if ~parse_ok
        warning('cpr_write_unit_spikes:parseError', ...
            'Could not parse unit_ID: "%s". Expected format: ch###_<tag>_unit#_range#. Skipping.', ...
            unit_id);
        continue
    end

    % --- Validate struct path -------------------------------------------
    % The spike data lives at phy.brain.CPR.spks.<channel>.<unit_field>{iCycle}
    if ~isfield(spks, channel)
        warning('cpr_write_unit_spikes:missingChannel', ...
            'Channel field "%s" not found in spks struct. Skipping %s.', ...
            channel, unit_id);
        continue
    end
    ch_struct = spks.(channel);

    if ~isfield(ch_struct, unit_field)
        warning('cpr_write_unit_spikes:missingUnit', ...
            'Unit field "%s" not found in spks.%s. Skipping %s.', ...
            unit_field, channel, unit_id);
        continue
    end
    unit_cell = ch_struct.(unit_field);   % {1 x nTotalCycles} spike time cells

    % --- Cycle mask for this range --------------------------------------
    cyc_ids = tbl.cyc_id{iRow};                   % [nCycles x 1] double
    cyc_ids = sort(cyc_ids(~isnan(cyc_ids)));      % remove NaN, sort ascending

    if isempty(cyc_ids)
        warning('cpr_write_unit_spikes:emptyCycMask', ...
            'cyc_id is empty for unit %s. Skipping.', unit_id);
        continue
    end

    % --- Collect spikes across all valid cycles for this range ----------
    all_cycle_ids   = [];   % [nSpikesTotal x 1] — grows with each cycle
    all_spike_times = [];   % [nSpikesTotal x 1] — microseconds

    for iCyc = 1:numel(cyc_ids)
        c = cyc_ids(iCyc);   % 1-based cycle number

        % Guard against cycle index exceeding the cell array length
        if c > numel(unit_cell)
            warning('cpr_write_unit_spikes:cycleOutOfRange', ...
                'Unit %s: cycle %d exceeds cell array length (%d). Skipping cycle.', ...
                unit_id, c, numel(unit_cell));
            continue
        end

        spk_times = unit_cell{c};   % [nSpikes x 1] double, microseconds

        % A cycle with zero spikes contributes no rows (not an error)
        if isempty(spk_times)
            continue
        end

        n = numel(spk_times);
        all_cycle_ids   = [all_cycle_ids;   repmat(c, n, 1)];  %#ok<AGROW>
        all_spike_times = [all_spike_times; spk_times(:)    ];  %#ok<AGROW>
    end

    % Skip if no spikes found across all cycles (unusual but possible)
    if isempty(all_spike_times)
        warning('cpr_write_unit_spikes:noSpikes', ...
            'No spikes found for %s across %d cycles. No file written.', ...
            unit_id, numel(cyc_ids));
        continue
    end

    % --- Write CSV -------------------------------------------------------
    out_name = sprintf('%s_%s_%s_%s_spikes.csv', rec, channel, unit_field, range_str);
    out_path = fullfile(spikes_dir, out_name);

    write_spike_csv(out_path, all_cycle_ids, all_spike_times, ...
        unit_id, channel, unit_field, range_str, cyc_ids, rec, tbl(iRow,:));

    n_written = n_written + 1;
    fprintf('[cpr_write_unit_spikes] %s  ->  %s  (%d spikes, %d cycles)\n', ...
        unit_id, out_name, numel(all_spike_times), numel(cyc_ids));
end

fprintf('[cpr_write_unit_spikes] Done. %d file(s) written to: %s\n', ...
    n_written, spikes_dir);

end % main function


% =========================================================================
%  LOCAL HELPER FUNCTIONS
% =========================================================================

function [channel, unit_field, range_str, ok] = parse_unit_id(unit_id)
% PARSE_UNIT_ID  Extract channel, unit field, and range from a unit_ID string.
%
%   Expected format: "<channel>_<unit_field>_<range>"
%       e.g. "ch001_neg_unit1_range1"  ->  ch001_neg  /  unit1  /  range1
%
%   The regex captures the longest possible channel name before _unit#_range#,
%   accommodating channels with multiple underscore-separated components.
%
%   OUTPUT ok = false if the format does not match the expected pattern.

tokens = regexp(unit_id, '^(.+)_(unit\d+)_(range\d+)$', 'tokens', 'once');

if isempty(tokens)
    channel    = '';
    unit_field = '';
    range_str  = '';
    ok         = false;
else
    channel    = tokens{1};    % e.g. 'ch001_neg'
    unit_field = tokens{2};    % e.g. 'unit1'
    range_str  = tokens{3};    % e.g. 'range1'
    ok         = true;
end

end % parse_unit_id


% -------------------------------------------------------------------------

function write_spike_csv(out_path, cycle_ids, spike_times, ...
    unit_id, channel, unit_field, range_str, cyc_ids, rec, row)
% WRITE_SPIKE_CSV  Write comment header block + spike data to a CSV file.
%
%   The file opens with '#'-prefixed comment lines containing unit metadata.
%   These are skipped automatically by pandas, R, and MATLAB readtable.
%   The data section follows as a standard two-column CSV.
%
%   INPUTS
%       out_path        Full path to output file
%       cycle_ids       [nSpikes x 1] cycle number for each spike
%       spike_times     [nSpikes x 1] spike times in microseconds
%       unit_id, channel, unit_field, range_str — parsed identifiers
%       cyc_ids         [nCycles x 1] valid cycle numbers for this range
%       rec             Recording ID string
%       row             Single-row table slice from the include table

gen_date = datestr(now, 'yyyy-mm-dd');

fid = fopen(out_path, 'w');
if fid == -1
    error('cpr_write_unit_spikes:fileOpenFailed', ...
        'Could not open file for writing: %s', out_path);
end

% --- Comment header (human-readable; skipped by CSV parsers) -------------
fprintf(fid, '# CPR spike times — %s\n', rec);
fprintf(fid, '# Generated   : %s\n', gen_date);
fprintf(fid, '# AI source   : Claude (Anthropic) — https://claude.ai\n');
fprintf(fid, '#\n');
fprintf(fid, '# unit_ID     : %s\n', unit_id);
fprintf(fid, '# channel     : %s\n', channel);
fprintf(fid, '# unit_field  : %s\n', unit_field);
fprintf(fid, '# range       : %s\n', range_str);
fprintf(fid, '#   A "range" is the span of cycles over which this unit was\n');
fprintf(fid, '#   tracked as a stable spike cluster. The same channel+unit\n');
fprintf(fid, '#   across different ranges is the same sorted cluster.\n');
fprintf(fid, '#   Cycle spans: %d to %d  (%d cycles total)\n', ...
    min(cyc_ids), max(cyc_ids), numel(cyc_ids));
fprintf(fid, '#\n');
fprintf(fid, '# Inclusion statistics (stimulus-driven test):\n');
fprintf(fid, '#   mean_FR           : %.4f Hz\n',    row.mean_FR);
fprintf(fid, '#   signrank_p        : %.4e\n',        row.signrank_p);
fprintf(fid, '#   signrank_z        : %.4f\n',        row.signrank_z);
fprintf(fid, '#   thresh_latency_ms : %.1f ms\n',     row.thresh_latency_ms);
fprintf(fid, '#\n');
fprintf(fid, '# n_spikes    : %d\n', numel(spike_times));
fprintf(fid, '# n_cycles    : %d\n', numel(unique(cycle_ids)));
fprintf(fid, '#\n');
fprintf(fid, '# COLUMNS\n');
fprintf(fid, '#   cycle_id       Cycle number (1-based; matches *_cycle<NNN>_frames.csv)\n');
fprintf(fid, '#   spike_time_us  Spike time in microseconds, relative to cycle onset\n');
fprintf(fid, '#                  Divide by 1e6 for seconds (aligns with time_s in frame CSV)\n');
fprintf(fid, '#\n');
fprintf(fid, '# FORMAT NOTE: Long format — one row per spike.\n');
fprintf(fid, '#   Cycles in the valid range with zero spikes produce no rows.\n');
fprintf(fid, '#   Load with:\n');
fprintf(fid, '#     Python : pd.read_csv(f, comment=''#'')\n');
fprintf(fid, '#     MATLAB : readtable(f, ''CommentStyle'', ''#'')\n');
fprintf(fid, '#     R      : read.csv(f, comment.char=''#'')\n');
fprintf(fid, '#\n');

% --- CSV column header and data ------------------------------------------
fprintf(fid, 'cycle_id,spike_time_us\n');
fprintf(fid, '%d,%d\n', [cycle_ids(:), round(spike_times(:))]');

fclose(fid);

end % write_spike_csv

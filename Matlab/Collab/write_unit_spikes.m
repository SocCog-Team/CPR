function write_unit_spikes(phy, output_dir)
% WRITE_UNIT_SPIKES  Export spike times for significantly stimulus-driven units.
%
%   WRITE_UNIT_SPIKES(PHY, OUTPUT_DIR) reads the spike inclusion table and
%   the spike time structs, then writes one long-format CSV per physical unit.
%
% -------------------------------------------------------------------------
%   INCLUSION CRITERIA
%       Only units where inclusion_flag == true are exported.
%       Rows sharing the same channel and unit field but different range
%       suffixes (e.g. unit1_range1, unit1_range2) represent the same
%       physical neuron recorded across different cycle spans. These are
%       merged into a single output file.
%
%   CYCLE SELECTION
%       For each unit row, phy.brain.CPR.spks.include.cyc_id specifies the
%       exact cycle numbers where the unit was present. Only those cycles
%       are exported. Cycles not listed in cyc_id are silently skipped.
%
% -------------------------------------------------------------------------
%   OUTPUT FORMAT  (long-format CSV)
%
%       One file per physical unit:
%           <rec_num>_<channel>_<unit_field>_spikes.csv
%
%       Columns:
%           cycle_id        Cycle number (1-based, matches cycle CSV filenames)
%           spike_time_us   Spike time in microseconds, relative to cycle onset
%
%       Each row is one spike. Cycles present in cyc_id but containing zero
%       spikes produce no rows (their absence is informative).
%
%       Example row:   23, 145230
%
% -------------------------------------------------------------------------
%   STRUCT ACCESS PATTERN
%       unit_ID "ch001_neg_unit1_range1" maps to:
%           phy.brain.CPR.spks.ch001_neg.unit1{iCycle}
%       The range suffix is used only to identify which cyc_id rows to use.
%
% -------------------------------------------------------------------------
%   REQUIREMENTS
%       MATLAB R2019b or later.
%
% -------------------------------------------------------------------------
%   EXAMPLE
%       write_unit_spikes(phy, fullfile('data', phy.exp.rec_num))
%
%   See also: write_cycle_csv, write_experiment_headers

% =========================================================================
% 0.  SETUP
% =========================================================================

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

tbl   = phy.brain.CPR.spks.include;        % 301×8 table
spks  = phy.brain.CPR.spks;                % struct with channel substructs
rec   = phy.exp.rec_num;                   % e.g. 'rec068'

% Keep only significantly stimulus-driven units
included_rows = find(tbl.inclusion_flag);
fprintf('[write_unit_spikes] %d / %d units pass inclusion criterion.\n', ...
    numel(included_rows), height(tbl));

% =========================================================================
% 1.  PARSE unit_IDs AND GROUP RANGES INTO PHYSICAL UNITS
%
%   unit_ID format:  "<channel>_<unit_field>_range<N>"
%   e.g.  "ch001_neg_unit1_range1"  ->  channel="ch001_neg", unit="unit1"
%
%   Multiple rows with the same (channel, unit_field) but different range
%   numbers are merged: their cyc_id lists are concatenated and deduplicated.
% =========================================================================

% Map from base_id string -> struct with channel, unit_field, and cycle list
unit_map = containers.Map();

for iRow = included_rows(:)'

    raw_id = char(tbl.unit_ID(iRow));       % e.g. 'ch001_neg_unit1_range1'
    cyc_id = tbl.cyc_id{iRow};             % [nCycles x 1] double

    % Strip the trailing _range<N> to get the base physical unit ID
    base_id = regexprep(raw_id, '_range\d+$', '');  % 'ch001_neg_unit1'

    % Parse channel and unit_field from base_id
    % Convention: last token is the unit field, everything before is channel
    tokens = strsplit(base_id, '_');        % {'ch001', 'neg', 'unit1'}
    unit_field   = tokens{end};             % 'unit1'
    channel_name = strjoin(tokens(1:end-1), '_');  % 'ch001_neg'

    % Accumulate cycle IDs across ranges for the same physical unit
    if isKey(unit_map, base_id)
        entry         = unit_map(base_id);
        entry.cyc_ids = union(entry.cyc_ids, cyc_id(~isnan(cyc_id)));
        unit_map(base_id) = entry;
    else
        entry.channel    = channel_name;
        entry.unit_field = unit_field;
        entry.base_id    = base_id;
        entry.cyc_ids    = cyc_id(~isnan(cyc_id));
        unit_map(base_id) = entry;
    end
end

% =========================================================================
% 2.  EXPORT ONE CSV PER PHYSICAL UNIT
% =========================================================================

unit_keys  = keys(unit_map);
n_units    = numel(unit_keys);
n_exported = 0;

for iUnit = 1:n_units

    entry      = unit_map(unit_keys{iUnit});
    ch         = entry.channel;        % 'ch001_neg'
    uf         = entry.unit_field;     % 'unit1'
    cyc_ids    = sort(entry.cyc_ids);  % sorted cycle numbers for this unit

    % Validate struct path before accessing
    if ~isfield(spks, ch)
        warning('write_unit_spikes:missingChannel', ...
            'Channel field "%s" not found in spks struct. Skipping.', ch);
        continue
    end
    ch_struct = spks.(ch);

    if ~isfield(ch_struct, uf)
        warning('write_unit_spikes:missingUnit', ...
            'Unit field "%s" not found in spks.%s. Skipping.', uf, ch);
        continue
    end
    unit_cell = ch_struct.(uf);   % {1 x nTotalCycles} cell of spike time vectors

    % ------------------------------------------------------------------
    % 2a. Collect spikes across all valid cycles
    % ------------------------------------------------------------------
    all_cycle_ids   = [];   % will grow: [nSpikesTotal x 1]
    all_spike_times = [];   % will grow: [nSpikesTotal x 1]

    for iCyc = 1:numel(cyc_ids)
        c = cyc_ids(iCyc);   % 1-based cycle number

        % Guard: cycle index must exist in the cell array
        if c > numel(unit_cell)
            warning('write_unit_spikes:cycleOutOfRange', ...
                'Unit %s: cycle %d exceeds cell array length (%d). Skipping.', ...
                entry.base_id, c, numel(unit_cell));
            continue
        end

        spk_times = unit_cell{c};   % [nSpikes x 1] double, microseconds

        % Skip empty cycles (no spikes recorded in this cycle)
        if isempty(spk_times)
            continue
        end

        n = numel(spk_times);
        all_cycle_ids   = [all_cycle_ids;   repmat(c, n, 1)];  %#ok<AGROW>
        all_spike_times = [all_spike_times; spk_times(:)    ];  %#ok<AGROW>
    end

    % Skip units with no spikes at all (edge case)
    if isempty(all_spike_times)
        warning('write_unit_spikes:noSpikes', ...
            'Unit %s has no spikes in any valid cycle. No file written.', ...
            entry.base_id);
        continue
    end

    % ------------------------------------------------------------------
    % 2b. Write CSV
    % ------------------------------------------------------------------
    out_name = sprintf('%s_%s_%s_spikes.csv', rec, ch, uf);
    out_path = fullfile(output_dir, out_name);

    write_spike_csv(out_path, all_cycle_ids, all_spike_times, entry, rec);

    n_exported = n_exported + 1;
    fprintf('[write_unit_spikes] %s  ->  %s  (%d spikes, %d cycles)\n', ...
        entry.base_id, out_name, numel(all_spike_times), numel(cyc_ids));
end

fprintf('[write_unit_spikes] Done. %d unit file(s) written to: %s\n', ...
    n_exported, output_dir);

end % main function


% =========================================================================
% LOCAL HELPER
% =========================================================================

function write_spike_csv(out_path, cycle_ids, spike_times, entry, rec)
% WRITE_SPIKE_CSV  Write header comment block and data to a CSV file.
%
%   The file begins with comment lines (prefixed '#') describing the
%   content and column schema, followed by a standard CSV header row and
%   the data. Comment lines are ignored by most CSV parsers (e.g. pandas
%   read_csv with comment='#', MATLAB readtable with CommentStyle='#').

fid = fopen(out_path, 'w');
if fid == -1
    error('write_unit_spikes:fileOpenFailed', ...
        'Could not open file for writing: %s', out_path);
end

% --- Comment header (human-readable metadata, ignored by CSV parsers) ---
fprintf(fid, '# CPR spike times — %s\n', rec);
fprintf(fid, '# unit        : %s_%s\n', entry.channel, entry.unit_field);
fprintf(fid, '# n_cycles    : %d\n',    numel(unique(cycle_ids)));
fprintf(fid, '# n_spikes    : %d\n',    numel(spike_times));
fprintf(fid, '#\n');
fprintf(fid, '# COLUMNS\n');
fprintf(fid, '#   cycle_id       Cycle number (1-based, matches *_cycle<NN>.csv filenames)\n');
fprintf(fid, '#   spike_time_us  Spike time in microseconds, relative to cycle onset\n');
fprintf(fid, '#\n');
fprintf(fid, '# NOTE: Long format. Each row is one spike.\n');
fprintf(fid, '#       Cycles present in cyc_id but with zero spikes have no rows.\n');
fprintf(fid, '#       To align with the frame-wise CSV, divide spike_time_us by 1e6\n');
fprintf(fid, '#       to obtain spike_time_s, comparable to the time_s column.\n');
fprintf(fid, '#\n');

% --- CSV column header and data -----------------------------------------
fprintf(fid, 'cycle_id,spike_time_us\n');
fprintf(fid, '%d,%d\n', [cycle_ids(:), round(spike_times(:))]');

fclose(fid);

end % write_spike_csv

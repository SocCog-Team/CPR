%%% Reproduce human performance as replay for dyadic experiments

% Import frame-wise behavior
load('/Users/fschneider/Desktop/pilot_tab/20250812_fih_psy4.mat')

% init
rdp_dir = nan([7200 length(out.raw.js_dir)]);
rdp_coh = nan([7200 length(out.raw.js_dir)]);
js_dir = nan([7200 length(out.raw.js_dir)]);
js_tlt = nan([7200 length(out.raw.js_dir)]);

for iCyc = 1:length(out.raw.js_dir)
 rdp_dir(1:length(out.raw.js_dir{iCyc}),iCyc) = out.raw.rdp_dir{iCyc};
 rdp_dir(:,iCyc) = replace_nan_zero(rdp_dir(:,iCyc));

 js_dir(1:length(out.raw.js_dir{iCyc}),iCyc) = out.raw.js_dir{iCyc};
 js_dir(:,iCyc) = replace_nan_zero(js_dir(:,iCyc));

 js_tlt(1:length(out.raw.js_dir{iCyc}),iCyc) = out.raw.js_ecc{iCyc};
 js_tlt(:,iCyc) = replace_nan_zero(js_tlt(:,iCyc));
 
 coh_vec = ((iCyc * 6)-5):1:(iCyc * 6);
 tmp_coh = [];
 
for i = coh_vec
     tmp_coh = [tmp_coh; repmat(out.coherence(i),1200,1)];
end

rdp_coh(:,iCyc) = tmp_coh;
end


% Export to summary file
% Produce scipt to export relevant data for each trial in .txt format


function out = replace_nan_zero(in)
v = in;

% Find the first valid entry (non-NaN, non-zero)
first_valid_idx = find(~isnan(v) & v ~= 0, 1, 'first');
first_valid_value = v(first_valid_idx);

% Replace leading NaNs with first valid value
v(1:first_valid_idx-1) = first_valid_value;

% Find the last valid entry (non-NaN, non-zero)
last_valid_idx = find(~isnan(v) & v ~= 0, 1, 'last');
last_valid_value = v(last_valid_idx);

% Replace trailing zeros with last valid value
v(last_valid_idx+1:end) = last_valid_value;
out = v;
end
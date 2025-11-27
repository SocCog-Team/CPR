
addpath /Users/fschneider/Documents/GitLab/matlab4mworks/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/Helper_functions
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/

load(['/Users/fschneider/ownCloud/Documents/Conferences/SfN_2025/data/summary_' rec_lst{1} '.mat'])
load(['/Users/fschneider/ownCloud/Documents/Conferences/SfN_2025/data/state_responses_' rec_lst{1} '.mat'])

state.include
plot_spkjoy_correlation(state,'ch015_neg_unit1','/Users/fschneider/Desktop/')

function plot_spkjoy_correlation(state,unit_id,dest_dir)
snr                 = unique(state.rdp_coh);
col                 = cool(length(snr));
FR                  = state.spk_n.(unit_id) ./ state.dur_s;
bin_width           = 30;
incl                = state.include.(unit_id);
dyad_idx            = contains(state.task,'dyad');
solo_idx            = contains(state.task,'solo');

for iState = 1:length(state.dur_s)
    js_err              = rad2deg(circ_dist(deg2rad(state.js_monk_dir{iState}),deg2rad(state.rdp_dir(iState)))); % Get circular distance to RDP direction
    state.js_acc{iState}= abs(1 - abs(js_err / 180));                               % Calculate accuracy
end

acc                 = cellfun(@mean,state.js_acc);
tlt                 = cellfun(@mean,state.js_monk_tlt);

[PD,~,VS] = preferredDirection(state.rdp_dir(incl), FR(incl));
roi                 = mod([PD-(bin_width/2) PD+(bin_width/2)],360);

if roi(1) > roi(2)
    PD_bin          = state.rdp_dir > roi(1) | state.rdp_dir < roi(2);
else
    PD_bin          = state.rdp_dir > roi(1) & state.rdp_dir < roi(2);
end

[p,~,stats] = ranksum(FR(PD_bin & solo_idx & incl), FR(PD_bin & ~dyad_idx & incl));  % Mann–Whitney–U-Test (nichtparametrisch)
% z = stats.zval;
[mean(FR(PD_bin & incl)) mean(FR(~PD_bin & incl))]

%%% PLOT %%%
figure; hold on
for iCoh = 1:length(snr)
    subplot(2,2,iCoh); hold on

    coh_idx = state.rdp_coh == snr(iCoh);
    sc = scatter(tlt(coh_idx & PD_bin & incl),FR(coh_idx & PD_bin & incl));
    sc.Marker = '.';
    sc.SizeData = 50;
    sc.MarkerEdgeColor = col(iCoh,:);
    sc.SizeData = 50;

    [pfit] = polyfit(tlt(coh_idx & PD_bin & incl),FR(coh_idx & PD_bin & incl),1);
    yfit = polyval(pfit,[0:.1:1]);
    plot([0:.1:1],yfit,'Color', col(iCoh,:),'LineWidth',2);
    [r,p] = corrcoef(tlt(coh_idx & PD_bin & incl),FR(coh_idx & PD_bin & incl));
    title(['r: ' num2str(r(2)) ' | p: ' num2str(p(2))])
    set(gca,'fontsize',20)
    xlabel('JS tilt [norm]')
    ylabel('FR [Hz]')
    axis square
end
% print([dest_dir '/corr_spk_joy_tlt_' char(unit_id)],'-dsvg')

figure; hold on
for iCoh = 1:length(snr)
    subplot(2,2,iCoh); hold on
    coh_idx =  state.rdp_coh == snr(iCoh);
    sc = scatter(acc(coh_idx & PD_bin & incl),FR(coh_idx & PD_bin & incl));
    sc.Marker = '.';
    sc.SizeData = 50;
    sc.MarkerEdgeColor = col(iCoh,:);
    sc.SizeData = 50;
    [pfit] = polyfit(acc(coh_idx & PD_bin & incl),FR(coh_idx & PD_bin & incl),1);
    yfit = polyval(pfit,[0:.1:1]);
    plot([0:.1:1],yfit,'Color', col(iCoh,:),'LineWidth',2);
    [r,p] = corrcoef(acc(coh_idx & PD_bin & incl),FR(coh_idx & PD_bin & incl));
    title(['r: ' num2str(r(2)) ' | p: ' num2str(p(2))])
    xlabel('JS accuracy [norm]')
    ylabel('FR [Hz]')
    set(gca,'fontsize',20)
    axis square
end
% print([dest_dir '/corr_spk_joy_acc_' char(unit_id)],'-dsvg')
end

function [prefDir, prefMag, vecStrength] = preferredDirection(angles, rates)
% preferredDirection computes the resultant vector and tuning strength from polar data.
%
%   [prefDir, prefMag, vecStrength] = preferredDirection(angles, rates)
%
%   INPUTS:
%       angles - vector of stimulus directions (in degrees or radians)
%       rates  - vector of corresponding firing rates
%
%   OUTPUTS:
%       prefDir     - preferred direction (same unit as input angle)
%       prefMag     - magnitude of the resultant vector
%       vecStrength - normalized vector magnitude (0–1), i.e. direction selectivity
%
%   Example:
%       angles = 0:45:315;
%       rates  = [5 8 12 9 4 3 2 6];
%       [dir, mag, vs] = preferredDirection(angles, rates)
%
%   See also: atan2, deg2rad, rad2deg

% Check input size
if numel(angles) ~= numel(rates)
    error('angles and rates must have the same length.');
end

% Detect whether angles are in degrees or radians
if max(abs(angles)) > 2*pi
    angRad = deg2rad(angles);
    useDegrees = true;
else
    angRad = angles;
    useDegrees = false;
end

% Compute resultant vector components
x = sum(rates .* cos(angRad));
y = sum(rates .* sin(angRad));

% Resultant vector magnitude
prefMag = sqrt(x^2 + y^2);

% Preferred direction (angle of resultant)
prefDir = atan2(y, x);

% Convert to degrees if needed
if useDegrees
    prefDir = rad2deg(prefDir);
    prefDir = mod(prefDir, 360);
end

% Normalized vector strength (0–1)
vecStrength = prefMag / sum(rates);
end

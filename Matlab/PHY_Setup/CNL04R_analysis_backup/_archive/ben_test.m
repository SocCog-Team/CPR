addpath('/Users/cnl/Desktop/CPR/PlexonMatlabOfflineFilesSDK/');
clear resolution minmax minmax_cnt

% cd /Users/fschneider/Desktop
cd /Users/cnl/Documents/DATA/Nilan/pl2/

% pl2_name = '/Volumes/T7_Shield/plexon/other/noise_1x_250gain.pl2';
pl2_name = '20250806_nil_CPR_block1_phy4_rec044_ann.pl2';
% pl2_name = '20250718_nil_CPR_block1_phy4_rec038_ann.pl2';
% pl2_name = '20250717_nil_CPR_block1_phy4_rec037_ann.pl2'; % amp_gain_x20
pl2         = PL2GetFileIndex(pl2_name);

for iChan = 1%64   
    clear ad abs_dff_sorted
    [ad]                = PL2Ad(pl2_name, iChan);
    abs_dff_sorted      = unique(sort(abs(diff(ad.Values))));
    resolution(iChan,:) = abs_dff_sorted(1:2);         
    minmax(iChan,:)     = [min(ad.Values) max(ad.Values)];
    minmax_cnt(iChan,:) = [sum(ad.Values == min(ad.Values)) sum(ad.Values == max(ad.Values))];
end

% save(['/Users/fschneider/Desktop/minmax_1x_250gain'],'minmax','-v7.3')
% save(['/Users/fschneider/Desktop/20250718_minmax_20x_250gain'],'minmax','-v7.3')
% save(['/Users/fschneider/Desktop/20250717_minmax_1x_250gain'],'minmax','-v7.3')
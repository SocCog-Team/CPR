addpath '/Volumes/cnl/Desktop/CPR/code'
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab

data_pth        ='/Volumes/cnl/Documents/MWorks/Data/';
cfg_pth         ='/Volumes/cnl/Desktop/CPR/code/felix_nhp_solo.cfg';
writeH5_flag    = true;
dest_dir_tbl    ='/Volumes/cnl/Desktop/CPR/mat/';
dest_dir_plot   ='/Volumes/cnl/Desktop/CPR/png/';

mwk2_lst = {'20240903_Dunkin_CPRsolo_block1_physio4.mwk2',...
    '20240905_Nilan_CPRsolo_block1_physio4.mwk2','20240905_Dunkin_CPRsolo_block1_physio4.mwk2'};

mwk2_lst = {'20240903_Dunkin_CPRsolo_block1_physio4.mwk2'}
for iFile = 1:length(mwk2_lst)
    fname = mwk2_lst{iFile};
    [d,t] = PHY_daily_summary_test(data_pth, fname, cfg_pth, writeH5_flag, dest_dir_tbl, dest_dir_plot);
end
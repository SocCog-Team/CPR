 % Setup
addpath /Users/cnl/Desktop/CPR/code
addpath /Users/cnl/Documents/GitLab/matLab4mworks
addpath /Users/cnl/Desktop/CPR/CircStat2012a
addpath /Users/cnl/Desktop/CPR/Violinplot-Matlab

data_pth        ='/Users/cnl/Documents/DATA/Nilan/';
cfg_pth         ='/Users/cnl/Desktop/CPR/code/felix_nhp_solo.cfg';
dest_dir_tbl    ='/Users/cnl/Desktop/CPR/mat/';
dest_dir_plot   ='/Users/cnl/Desktop/CPR/png/';

% Server
% addpath /Users/cnl/Desktop/CPR/MatLab4MWorks
% addpath /Users/cnl/Desktop/CPR/CPR/CircStat2012a
% data_pth        ='/Users/cnl/Documents/MWorks/Data/';
% cfg_pth         ='/Users/cnl/Desktop/CPR/code/felix_nhp_solo.cfg';
% writeH5_flag    = false;
% dest_dir_tbl    ='/Users/cnl/Desktop/CPR/mat/';
% dest_dir_plot   ='/Users/cnl/Desktop/CPR/png/';

mwk2_lst        = {'20250807_nil_CPR_block1_phy4_rec045_ann.mwk2'};
writeH5_flag    = true;

for iFile = 1:length(mwk2_lst)
    fname       = mwk2_lst{iFile};
    [d,t]       = PHY_daily_summary_test(data_pth, fname, cfg_pth, writeH5_flag, dest_dir_tbl, dest_dir_plot);
end
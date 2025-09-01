% addpath '/Volumes/cnl/Desktop/CPR/code'
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab

% data_pth        ='/Volumes/cnl/Documents/MWorks/Data/';
% cfg_pth         ='/Volumes/cnl/Desktop/CPR/code/felix_nhp_solo.cfg';
% writeH5_flag    = true;
% dest_dir_tbl    ='/Volumes/cnl/Desktop/CPR/mat/';
% dest_dir_plot   ='/Volumes/cnl/Desktop/CPR/png/';

% mwk2_lst = {'20240903_Dunkin_CPRsolo_block1_physio4.mwk2',...
%     '20240905_Nilan_CPRsolo_block1_physio4.mwk2','20240905_Dunkin_CPRsolo_block1_physio4.mwk2'};

data_pth        ='/Users/fschneider/Desktop/A01/mwk2/';
cfg_pth         ='/Users/fschneider/Documents/GitHub/CPR/Matlab/PHY_Setup/server_code/felix_nhp_solo.cfg';
writeH5_flag    = true;
dest_dir_tbl    ='/Users/fschneider/Desktop/A01/tbl/';
dest_dir_plot   ='/Users/fschneider/Desktop/A01/png/';

% mwk2_lst = {'20250312_nil_CPRsolo_block1_phy4.mwk2'};
cd(data_pth)
mwk2_lst = dir('*.mwk2');

%%
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/PHY_Setup/server_code/
for iFile = 4:length(mwk2_lst)
    disp(['File: ' num2str(iFile) '/' num2str(length(mwk2_lst))])
%     fname = mwk2_lst{iFile};
    fname = mwk2_lst(iFile).name;
    [d,t] = copy_PHY_daily_summary_test(data_pth, fname, cfg_pth, writeH5_flag, dest_dir_tbl, dest_dir_plot);
end
%%% DATA CONVERSION WOLF LAB %%%

close all
clear all

% Import preprocessed example data
load('/Users/fschneider/Desktop/summary_20251204_nil_CPR_block1_phy4_rec068_fxs.mat')

%% Export experimental cycle header

output_dir = '/Users/fschneider/Desktop/rec068/';
file_prefix = sprintf('%s_%s_%s_%s', ...
    phy.exp.date, phy.exp.monkey, phy.exp.rec_num, phy.exp.block);

write_experiment_headers(phy, file_prefix, output_dir)

for iCycle = 1:length(phy.stim.rdp_dir)
    write_cycle_csv(phy, iCycle, file_prefix, output_dir)
end

write_unit_spikes(phy, file_prefix, output_dir)
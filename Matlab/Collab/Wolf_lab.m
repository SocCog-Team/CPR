%%% DATA CONVERSION WOLF LAB %%%

close all
clear all
clc

addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/Collab

% Import preprocessed example data
load('/Users/fschneider/Desktop/summary_20251204_nil_CPR_block1_phy4_rec068_fxs.mat')
%% Export experimental cycle header

root = fullfile('export', phy.exp.rec_num);

cpr_write_session_header(phy, root);          % README + cycle JSONs

for iCycle = 1:size(phy.stim.cpr_cyle, 1)
    cpr_write_cycle_frames(phy, iCycle, root);
end

cpr_write_unit_spikes(phy, root);
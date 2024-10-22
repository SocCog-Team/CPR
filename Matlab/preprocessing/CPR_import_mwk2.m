function d = CPR_import_mwk2(fname, var_lst, write_file, cfg_pth)

% This function imports .mwk2 files into a data structure using the 
% MW_readFile function developed by Ralf Brockhausen.
%
% Input:  	.fname          String or cell containing multiple strings.
%                           Needs to match with file.
%          	.var_lst        Cell, Contains strings of variables/variable
%                           classes that are imported.
%           .write_file     Boolean, Indicates if structure is written to
%                           .mat file in same directory.
%
% Output:   .d              Structure, Data structure containing variable 
%                           names, time stamps and respective values. 
%
%
% Example:	d = CPR_import_mwk2('fname', [], false)
%
% Known bugs:
%
% Feature wish list:
%
% Version history
%   1.0     (fxs 2020-09-01) Initial version.
%   1.1     (fxs 2021-11-24) HDF5 data writing.

%% Check input

addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab

if nargin < 4
    cfg_pth	= '/Users/fschneider/Documents/GitHub/CPR/Matlab/CFG/felix_solo_vs_dyadic.cfg';
end

if nargin < 3 || isempty(write_file)
    write_file  = false;
end

if nargin < 2 || isempty(var_lst)
    var_lst     = {
        'ML_', ...
        'CTRL_', ...
        'RDP_', ...
        'INFO_', ...
        'TRIAL_', ...
        'IO_joystickDirection', ...
        'IO_joystickStrength',...
        'IO_fixation_flag',...
        'EYE_x_dva',... 
        'EYE_y_dva'};
end

if nargin < 1 || isempty(fname)
   error('File name not specified') 
end

%% Import data file

if iscell(fname)                                                            % If multiple files ...
    if isfile([fname{1}(1:12) '_merged.mat']) && write_file == false   
        d           = MW_readH5([fname{1}(1:34) '_merged.h5']);
    else
        d.time      = [];                                                  	% Otherwise, import & merge...
        d.event     = [];
        d.value     = [];
        
        for i = 1:size(fname,2)
            tmp    	= MW_readFile(fname{i}, 'include', var_lst, '~typeOutcomeCheck','dotPositions');
            
            d.time  = [d.time tmp.time];
            d.event = [d.event tmp.event];
            d.value = [d.value tmp.value];
        end
        
        if write_file
            disp('Save struct...')
            MW_writeH5(d, [fname{1}(1:end-5) 'CPR_merged.h5'], 'replace', 'privateCFG', cfg_pth)  
            disp('Done!')
        end
    end
else
    if isfile([fname(1:end-5) '.h5']) && write_file == false                 % If file available...
        d       = MW_readH5([fname(1:end-5) '.h5']);                        % ...load .h5 file
    else
        d     	= MW_readFile(fname, 'include', var_lst, '~typeOutcomeCheck','dotPositions'); % Import .mwk2 file
        
        if write_file
            disp('Save data structure...')
            MW_writeH5(d, [fname(1:end-5) '.h5'], 'replace', 'privateCFG', cfg_pth) % Save to .h5
            disp('Done!')
        end
    end
end
end

function d = CPR_import_mwk2(fname, var_lst, write_file)

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

%% Check input

addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab

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
        'IO_joystickStrength' };
end

if nargin < 1 || isempty(fname)
   error('File name not specified') 
end

%% Import data file
if iscell(fname)                                                            % If multiple files ...
    if isfile([fname{1}(1:12) '_merged.mat'])
        tmp         = load([fname{1}(1:34) '_merged.mat']);                 % Lod if merged file already exists
        d           = tmp.d;
    else
        d.time      = [];                                                  	% Otherwise, import & merge...
        d.event     = [];
        d.value     = [];
        
        for i = 1:size(fname,2)
            tmp    	= MW_readFile(fname{i}, 'include', var_import);
            
            d.time  = [d.time tmp.time];
            d.event = [d.event tmp.event];
            d.value = [d.value tmp.value];
        end
        
        if write_file
            disp('Save struct...')
            save([fname{1}(1:12) '_CPR_merged.mat'], 'd', '-v7.3')
            disp('Done!')
        end
    end
else
    if isfile([fname '.mat'])                                               % If .mat file available...
        tmp   	= load([fname '.mat']);                                  	% ...load .mat file
        d   	= tmp.d;
    else
        d     	= MW_readFile(fname, 'include', var_lst);                   % Import .mwk2 sesion file
        
        if write_file
            disp('Save struct...')
            save([fname '.mat'], 'd', '-v7.3')                              % Save as .mat file
            disp('Done!')
        end
    end
end
end
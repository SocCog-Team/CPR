function out = getTrialData(dat,idx1,idx2)

% Helper function to extract MWorks trial data.
%
% Input:        .dat    MWorks input
%               .idx1   Trial index
%               .idx2   Variable index
%
% Outout:       .out    Corresponding data values
%
% Example:      out = getTrialData(d.value, trlIdx, varIdx)
%
% Felix Schneider, CNL, 11/2020

c = dat(idx1 & idx2);

if isempty(c) || sum(idx1 & idx2) == 0
    out	= nan;
else  
    if iscell(c)
        for i = 1:length(c)
            if isa(c{i},'int64')
                c{i} = double(c{i});
            end
        end
        
        if ischar(c{1})
        else
            c = cell2mat(c);
        end
    end
    out = c;
end
end
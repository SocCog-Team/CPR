function out = CPR_write_txt(STIM, pth)

% Check input
if nargin < 2
    pth             = [pwd '/'];
end

% Update text files with stimulus params
var_names           = {'RDP_direction',...
                        'RDP_coherence',...
                        'RDP_coherence_ts',...
                        'feedback_ts'};

for iVar = 1:size(var_names,2)
    clear vec fid
    fid             = fopen([pth var_names{iVar} '.txt'],'w');
    
    if iVar == 1
        vec         = sprintf('%d\n', STIM.RDP_direction);                  % Convert to tab-spaced vector
        vec         = vec(1:end-1);                                         % Remove final tab
    elseif iVar == 2
        vec         = sprintf('%f\n', STIM.RDP_coherence);
        vec         = vec(1:end-1); 
    elseif iVar == 3
        vec         = sprintf('%d\n', STIM.RDP_coherence_ts);
        vec         = vec(1:end-1);
    elseif iVar == 4
        vec         = sprintf('%d\n', STIM.feedback_ts);
        vec         = vec(1:end-1); 
    end
    
    fprintf(fid, vec);                                                      % Write to file
    fclose(fid);                                                            % Close file
end

% Update counter
fid                 = fopen([pth '/counter.txt'],'w');  
fprintf(fid, num2str(STIM.trial));
fclose(fid);

out                 = true;
end
function experimentData = MW_readH5(filename, varargin)
% function experimentData = MW_readH5(filename, varargin)
%
% Possible optional parameter:
%   'include', {list of keywords which should be read}
%       If you do not enter this parameter the function use the standard
%       includelist. In the includelist you can use the first letters, f.w.
%       'EYE_' for all the eyeparameter.
%
% Examples:
%   MW_readH5('bew-MPH5-phu-012.h5', 'include', {'TRIAL_', 'STIM_', 'IO_'})
%
% Known bugs:
%
% Feature wish list:
%
% Version history
%   1.0     (rab 2021-10-28) Initial version.
%   2.0     (rab 2023-06-07) New routine with high-level MatLab h5 calls.
%   2.0.1   (rab 2023-06-16) Convert strings to characterarrays


includeEvents = [];
skipArgument = false;

% go through the inputparameter
for (argCX = 1:size(varargin,2))    
    if ~skipArgument
        if (ischar(varargin{argCX}))
            %fprintf('     %s\n', varargin{i});
            switch varargin{argCX}
                case 'dotPositions'
                    dotPositions = true;
                case 'exclude'
                    dType = varargin{argCX+1}; % BENNENUNG? RABRABRAB
                    skipArgument = true;
                    excludeList = true;
                case 'include'
                    includeEvents = varargin{argCX+1};
                    skipArgument = true;
                    includeList = true;           
                case 'debugLevel'
                    debugLevel = varargin{argCX+1};
                    skipArgument = true; 
                otherwise
                    disp('WARNING: Unknown input argument *', varargin{argCX}, '*!');
            end
        elseif(isnumeric(varargin{argCX}))
            disp('WARNING: Unknown input argument *', varargin{argCX}, '*!');
        end
    else
        skipArgument = false;
    end
end


%filename = 'bew-MPH5-phu-012.h5'; % <- filename
%project = h5readatt(filename, '/', 'project');
%param2read  = {'/value', '/time'};
file_info = h5info(filename, '/value');

allEvents = {file_info.Datasets.Name};

if ~isempty(includeEvents)
        
    selectedEvents = 0;
    
    for inCX = 1:length(includeEvents)
        selectedEvents = selectedEvents | contains(allEvents, includeEvents{inCX});
    end
    events2read = allEvents(selectedEvents);
else
    events2read = {file_info.Datasets.Name};
end


sizeList = [];

% Um das Einlesen zu beschleunigen, ist es wichtig die ""großen" Brocken
% erst am Ende laden.
for fieldCX = 1:length(events2read)
    tempSize = file_info.Datasets(matches(allEvents, events2read{fieldCX})).Dataspace.Size;
    sizeList(fieldCX) = tempSize(length(tempSize));
end
[~, sizeIndex] = sort(sizeList);
events2read = events2read(sizeIndex);




tic
fprintf('MW_readH5: loading data...\n');



%% New "function"
data.value = {};
data.event = [];
data.time = [];


for fieldCX = 1:length(events2read)
    %events2read{fieldCX}
    readData.value = h5read(filename, ['/value/' events2read{fieldCX}]);
    readData.time = h5read(filename, ['/time/' events2read{fieldCX}]);
    
    if ( iscell(readData.value) )
        data.value = vertcat(data.value, readData.value);
    else
        if ~isnumeric(readData.value(1)) ...
            && isstr(readData.value{1})
                readData.value = char(readData.value);
                data.value =  vertcat(data.value, cellstr(readData.value));
            %end
        else
            data.value = vertcat(data.value, num2cell(readData.value, 2));
        end
    end
    data.time = vertcat(data.time, readData.time); %event_time.(events2read{fieldCX}));
    data.event = vertcat(data.event, repelem(string(events2read{fieldCX}), length(readData.time))');
     
end
data.event = categorical(data.event);
% don't forget to sort!

[~, order] = sort([data.time],'ascend');
data.time = data.time(order)';
data.event = data.event(order)';
data.value = data.value(order)';

experimentData = data;



%% Old
% 
% 
% 
% % read data into structs
% [event_value, event_time] = h5_extract(filename, '/', param2read, events2read);
% 
% experimentData.value = {};
% experimentData.event = [];
% experimentData.time = [];
% fList = fields(event_value);
% sizeList = [];
% 
% % Um das Einlesen zu beschleunigen, ist es wichtig die ""großen" Brocken
% % erst am Ende laden.
% for fieldCX = 1:length(fList)
%     fieldCX;
%     fList{fieldCX};
%     tempSize = file_info.Datasets(matches(allEvents, fList{fieldCX})).Dataspace.Size;
%     sizeList(fieldCX) = tempSize(length(tempSize));
%     %fprintf("%s\n", fList{fieldCX});
% end
% [~, sizeIndex] = sort(sizeList);
% fList = fList(sizeIndex);
% 
% 
% for fieldCX = 1:length(fList)
% 
%     fieldCX;
%     fList{fieldCX}
%     event_value.(fList{fieldCX})';
%     %cat(1, event_value.STIM_RDPblack_type, event_value.STIM_background_type , num2cell(event_value.STIM_RDPblack_coherence))
%     if ( iscell(event_value.(fList{fieldCX})) )
%         experimentData.value = horzcat(experimentData.value, event_value.(fList{fieldCX})');
%     else
%         experimentData.value = horzcat(experimentData.value, num2cell(event_value.(fList{fieldCX})));
%     end
%     experimentData.time = horzcat(experimentData.time, event_time.(fList{fieldCX})); % war '
%     experimentData.event = horzcat(experimentData.event, repelem(string(fList{fieldCX}),length(event_time.(fList{fieldCX}))));
%     %experimentData.value;
% end
% %dont forget to sort
% 
% experimentData.event = categorical(experimentData.event);
% toc
toc
end

% %% BIS HIIIIIIER
% 
% 
% times2probe = [2 * 10.^9, 3 * 10.^9];
% events2probe = events2read;
% probed_data = probe_time(event_value, event_time, events2probe, times2probe);
% 
% function probed_data = probe_time(event_value, event_time, events2probe, times2probe)
% % extracts parameters of experiment (events2probe) at speficified time points (times2probe)
% %   Input parameters:
% 
% %   event_value:        struct containing values of events
% %   event_time:         struct containing timing of events
% %   events2probe:       events to probe
% %   times2probe:        time points to probe, if unassigned, parameters at all
% %                       time points from event_time are probed
% 
% if nargin < 4
%     full_time_list = struct2cell(event_time);
%     times2probe = unique(cell2mat(full_time_list));
% end
% 
% probed_data = struct();
% for t = 1:length(times2probe)
%     zeit = times2probe(t);
%     for param = 1 : length(events2probe)
%         event_to_read = events2probe{param};
%         if sum(isspace(event_to_read)) > 0
%             event_to_read(isspace(event_to_read)) = '_';
%         end
%         times = event_time.(event_to_read);
%         position = find(zeit >= times, 1, 'last');
%         wert = event_value.(event_to_read)(position);
%         if length(times2probe) < 5
%             disp([zeit, {event_to_read}, wert(:)']);
%         end
%         probed_data.time(t) = zeit;
%         probed_data.(event_to_read){t} = wert;
%     end
% end
% end

%h5disp('20220428_jus_CPRsolo_block1_psycho4_fxs.h5','/')
%MW_readH5('20220428_jus_CPRsolo_block1_psycho4_fxs.h5', 'include', {'TRIAL_', 'STIM_', 'IO_','EYE_', 'INFO_'})
%hurz2 = h5info('20220428_jus_CPRsolo_block1_psycho4_fxs.h5')
function data = CPR_data_correction(data, var1, var2)

if (isstring(var1) || ischar(var1)) == false && (isstring(var2) || ischar(var2)) == false
    error('Variable names must be specified as strings or characters')
end

% Build index
idx.var1        = data.event == var1;
idx.var2        = data.event == var2;

% Extract timestamps and data
d.var1_ts   	= data.time(idx.var1);
d.var2_ts   	= data.time(idx.var2);
d.var1_val   	= data.value(idx.var1);
d.var2_val   	= data.value(idx.var2);

fprintf(['Correcting data entries for variables:\n' var1 '\n' var2 '\n']) 

% Find differences between sample timestamps
[val, pos]          = setdiff(d.var1_ts,d.var2_ts);

% Append missing data to input structure
for iPos = 1:length(pos)  
    if pos(iPos) > length(d.var2_ts)
        data.time	= [data.time d.var2_ts(end)];                    
        data.value  = [data.value d.var2_val(end)];
        data.event  = [data.event var2];
    else
        data.time	= [data.time d.var2_ts(pos(iPos))];                    
        data.value  = [data.value d.var2_val(pos(iPos))];
        data.event  = [data.event var2];
    end
end

% Check other way around
[val, pos]          = setdiff(d.var2_ts,d.var1_ts);

% Append missing data to input structure
for iPos = 1:length(pos)
    if pos(iPos) > length(d.var1_ts)
        data.time	= [data.time d.var1_ts(end)];
        data.value  = [data.value d.var1_val(end)];
        data.event  = [data.event var1];
    else
        data.time   = [data.time d.var1_ts(pos(iPos))];
        data.value  = [data.value d.var1_val(pos(iPos))];
        data.event  = [data.event var1];
    end
end

% Sort data
[data.time,i] = sort(data.time);
data.value = data.value(i);
data.event = data.event(i);

fprintf('Done!\n') 

end
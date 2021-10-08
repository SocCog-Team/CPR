function [aa,bb,ts] = CPR_adjust_vector_length(a,a_ts,b,b_ts)

aa                  = [];                                                   % Initiate output
bb                  = [];                                                   % Initiate output
ts                  = [];                                                   % Initiate output

if size(a,2) > size(b,2)
    vals            = setdiff(a_ts,b_ts);                                   % Set difference of two arrays
    pos             = [1 find(ismember(a_ts,vals))];                    	% Find positions
    
    n              	= diff([0 find(ismember(a_ts,vals)) length(b_ts)]);  	% Calculate segment length
    sub_val         = mat2cell(b,1,n);                                     	% Split into subvectors
    
    for i = 1:size(sub_val,2)-1
        sub_val{i}  = [sub_val{i} sub_val{i}(end)];                         % Copy missing value
    end
    
    aa              = a;
    bb              = cell2mat(sub_val);
    ts              = a_ts;
    
    bb(find(ismember(b_ts,setdiff(b_ts,a_ts)))) = [];                       % Delete other irregularities
    
else
    vals            = setdiff(b_ts,a_ts);                                   
    pos             = [1 find(ismember(b_ts,vals))];                    	
    
    n               = diff([0 find(ismember(b_ts,vals)) length(a_ts)]);   	
    sub_val         = mat2cell(a,1,n);                                     
    
    for i = 1:size(sub_val,2)-1
        sub_val{i}  = [sub_val{i} sub_val{i}(end)];                       	
    end
    
    aa              = cell2mat(sub_val);
    bb              = b;
    ts              = b_ts;
    
    aa(find(ismember(a_ts,setdiff(a_ts,b_ts)))) = [];
end

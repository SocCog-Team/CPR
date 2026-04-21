%%% Transform CPR data %%%
close all
clear all

% Source path
spth                    = '/Volumes/T7_Shield/CPR_psychophysics/';
folder                  = dir(spth);

for iFld = 1:length(folder)
    % Initialise
    t1                  = [];
    t2                  = [];
    fnames              = [];
    
    disp(['Folder: ' num2str(iFld) ' || ID: ' num2str(folder(iFld).name)])
    
    % Skip directories
    if iFld == 9 || iFld == 13 || iFld == 46 || iFld == 78 % weird scores - presumably restarted experiment?
        continue
    end
    
    if  strcmp(folder(iFld).name,'CaH') || ... % different experimental paradigm
            strcmp(folder(iFld).name,'DeM') || ...
            strcmp(folder(iFld).name,'DiH') || ...
            strcmp(folder(iFld).name,'GaA') || ...
            strcmp(folder(iFld).name,'MaM') || ...
            strcmp(folder(iFld).name,'MaS') || ...
            strcmp(folder(iFld).name,'NaK') || ...
            strcmp(folder(iFld).name,'SeH') || ...
            strcmp(folder(iFld).name,'ToS') || ...
            strcmp(folder(iFld).name,'Dyad18') || ...   % different coherence levels tested
            strcmp(folder(iFld).name,'Dyad19') || ...   % different coherence levels tested
            strcmp(folder(iFld).name,'Dyad36') || ...   % sessions aborted [eye signal]
            strcmp(folder(iFld).name,'Dyad55')
        continue
    end
    
    % Extract relevant file names
    if length(folder(iFld).name) == 3
        fnames          = getFiles([spth folder(iFld).name '/summary/'], 'agent');
    elseif length(folder(iFld).name) == 6
        fnames          = getFiles([spth folder(iFld).name '/summary/'], 'dyad');
    end
    
    if ~isempty(fnames)
        % 4 Files: 2 players x 2 Blocks
        if length(fnames) < 4
            error('Something is wrong...')
        end
        
        % Import preprocessed .mat data tables
        for iFile = 1:length(fnames)
            in          = load(fnames{iFile});
            
            if contains(fnames{iFile},'CPRagent')
                if contains(fnames{iFile},'agnt')
                    t1 	= [t1; in.t]; % Computer player
                else
                    t2 	= [t2; in.t]; % Human player
                end
            else
                % Files sorted alphabetically: [a a b b]
                if iFile <= length(fnames)/2
                    t1 	= [t1; in.t];
                else
                    t2 	= [t2; in.t];
                end
            end
        end
        
        % Arrange as 3D matrix
        dmat            = convert2matrix(t1,t2);
        
        % Crop NaNs - Reduce matrix to duration of shortest cycle
        out             = dmat(:,sum(isnan(dmat(1,:,:)),3)== 0,:);
        
        % Check if output if faulty
        check_output(out)
        
        % Write to file
        dest_dir    	= '/Users/fschneider/ownCloud/Wibral_lab/dyad_mat/';
        save([dest_dir 'dyad_' char(t1.ID(1)) '_' char(t2.ID(1))], 'out', '-v7.3')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = getFiles(pth, str)

cd(pth)
files                       = dir('*.mat');
c                           = 0;
out                         = [];

for iFile = 1:length(files)
    if contains(files(iFile).name, str)
        c                   = c+1;
        out{c}              = files(iFile).name;
    end
end
end

function out = convert2matrix(t1,t2)

out                         = [];

% For each block...
for iBlock = 1:2
    b_p1                    = t1(t1.block == iBlock,:);
    b_p2                    = t2(t2.block == iBlock,:);
    
%     % Convert cumulative- to individual target score
%     b_p1                  	= convertScore(b_p1);
%     b_p2                 	= convertScore(b_p2);
    
    % For each cycle...
    for iCycle = 1:b_p1.cyc_no(end)
        c_p1                = b_p1(b_p1.cyc_no == iCycle,:);
        c_p2                = b_p2(b_p2.cyc_no == iCycle,:);
        
        % Kick out last incomplete cycle
        if size(c_p1,1) < 30
            continue
        end
        
        % Initialise
        ts                  = [];
        rdp_dir             = [];
        vec_dir             = [];
        rdp_coh             = [];
        vec_len             = [];
        p1_dir              = [];
        p1_ecc              = [];
        p2_dir              = [];
        p2_ecc              = [];
        trg_ts            	= [];
        trg_score1        	= [];
        trg_score2        	= [];
        
        % Concatenate state-wise data to vector
        for iState = 1:size(c_p1,1)
            disp (['Block: ' num2str(iBlock) '/2 || Cycle: ' num2str(iCycle) '/' num2str(b_p1.cyc_no(end)) ' || State: ' num2str(iState) '/' num2str(size(c_p1,1))])
            
            nSamples        = length(c_p1.frme_ts{iState});
            ts              = [ts c_p1.frme_ts{iState}];
            rdp_dir         = [rdp_dir repmat(c_p1.rdp_dir(iState),1,nSamples)];
            rdp_coh         = [rdp_coh repmat(c_p1.rdp_coh(iState),1,nSamples)];
            
            [tmp_dir, tmp_len] = integrateDots(c_p1.rdp_dot,iState);
            
            if sum(isnan(tmp_dir)) == length(tmp_dir) || sum(isnan(tmp_len)) == length(tmp_len)
                % Add zeros if speed criterion was not met
                vec_dir    	= [vec_dir zeros(1,length(tmp_dir))];
                vec_len   	= [vec_len zeros(1,length(tmp_len))];
            else
                vec_dir    	= [vec_dir interpolateSamples(tmp_dir)];
                vec_len   	= [vec_len interpolateSamples(tmp_len)];
            end
            
            p1_dir          = [p1_dir c_p1.js_dir{iState}];
            p1_ecc          = [p1_ecc c_p1.js_ecc{iState}];
            p2_dir          = [p2_dir c_p2.js_dir{iState}];
            p2_ecc          = [p2_ecc c_p2.js_ecc{iState}];
            
            if c_p1.trg_shown(iState)
                trg_ts     	= [trg_ts c_p1.trg_ts{iState}];
                
                if sum(c_p1.trg_score{iState} > 1) > 0 || sum(c_p2.trg_score{iState} > 1) > 0
                    trg_score1	= [trg_score1 (c_p1.trg_acc{iState}.*c_p1.trg_ecc{iState}) .* c_p1.trg_hit{iState}];
                    trg_score2	= [trg_score2 (c_p2.trg_acc{iState}.*c_p2.trg_ecc{iState}) .* c_p2.trg_hit{iState}];
                else
                    trg_score1	= [trg_score1 c_p1.trg_score{iState}];
                    trg_score2	= [trg_score2 c_p2.trg_score{iState}];
                end
            end
        end
        
        % Add concatenated data to output matrix
        nan_mat             = nan(7,1e5);
        zer_mat             = zeros(5,1e5);
        mat                 = [nan_mat;zer_mat];
        nVec                = 1:length(ts);
        
        % Timestamps
        mat(1,nVec)         = ts;
        % Nominal stimulus direction
        mat(2,nVec)         = rdp_dir;
        % Actual stimulus direction based on dot movement
        mat(3,nVec)         = vec_dir;
        % Nominal stimulus coherence
        mat(4,nVec)         = rdp_coh;
        % Actual stimulus coherence based on dot movement [normalised]
        mat(5,nVec)         = vec_len;
        % Joystick direction: Player 1
        mat(6,nVec)         = interpolateSamples(mod(p1_dir,360)); % Interpolate NaNs in joystick samples
        % Joystick eccentricity: Player 1
        mat(7,nVec)         = interpolateSamples(p1_ecc);
        % Joystick direction: Player 2
        mat(8,nVec)         = interpolateSamples(mod(p2_dir,360));
        % Joystick eccentricity: Player 2
        mat(9,nVec)         = interpolateSamples(p2_ecc);     
        
        % Replace NaNs with zeros
        trg_score1(isnan(trg_score1)) = 0;
        trg_score2(isnan(trg_score2)) = 0;
        
        if sum(trg_score1 > 1) > 0 || sum(trg_score2 > 1) > 0
            error('shit!')
        end
        
        % Calculate running average [displayed on screen continuously]
        mavg_score1      	= runningAverage(trg_score1);
        mavg_score2      	= runningAverage(trg_score2);
        prev_trg_idx      	= 1;
        
        for iTrg = 1:length(trg_ts)
            % UNDERLYING ISSUE IN TARGET - STATE ASSIGNMENT -> Skip
            if sum(trg_ts(iTrg) > ts(end)) > 0
                warning(['target state assignment faulty: ' num2str([iBlock iCycle iTrg length(trg_ts)])])
                continue
            end
            
            trg_idx        	= find(ts >= trg_ts(iTrg),1,'first');
            % Timestamp of target appearance
            mat(10,trg_idx)	= trg_ts(iTrg);
            % Reward score player 1
            mat(11,trg_idx) = trg_score1(iTrg);
            % Reward score player 2
            mat(12,trg_idx)	= trg_score2(iTrg);
            
            if iTrg ~= length(trg_ts)
                pos      	= prev_trg_idx:trg_idx-1;
                n         	= trg_idx-prev_trg_idx;
            else
                pos         = prev_trg_idx:nVec(end);
                n       	= length(pos);
            end
            
            % Cumulative reward score: Player 1
            mat(13,pos)   	= repmat(mavg_score1(iTrg),1,n);
            % Cumulative reward score: Player 2
            mat(14,pos)   	= repmat(mavg_score2(iTrg),1,n);
            prev_trg_idx  	= trg_idx;
        end
        
        % Add to matrix to output
        out                 =  cat(3,out,mat);
    end
end
end

function out = interpolateSamples(in)

idx                         = find(isnan(in));
out                         = in;

for iSmpl = 1:length(idx)
    
    k                       = find(~isnan( out(idx(iSmpl):end) ),1,'first');
    n                       = idx(iSmpl) + (k-1);
    
    % First sample
    if idx(iSmpl) == 1
        out(idx(iSmpl))     = out(n);
        % Last sample
    elseif idx(iSmpl) == length(out)
        out(idx(iSmpl))     = length(out)-1;
        % Any other sample
    else
        out(idx(iSmpl))     = mean([out(idx(iSmpl)-1) out(n)]);
    end
end

% Check for NaNs
if sum(isnan(out)) > 0
    error('NaN!');
end
end

function out = runningAverage(in)
if sum(isnan(in))> 0
    error('NaN - not working')
end

for iTrg = 1:length(in)
    out(iTrg)               = sum(in(1:iTrg)) / iTrg;
end

out                         = [0 out];
end

function out = convertScore(in)

cum_score                 	= in.trg_score;
last_score              	= 0;

for iCell = 1:length(cum_score)
    for iTrg = 1:length(cum_score{iCell})
        new_score{iCell}(iTrg) = cum_score{iCell}(iTrg) - last_score;
        
        if new_score{iCell}(iTrg) < 0 || new_score{iCell}(iTrg) > 1
            disp([iCell iTrg])
            error('Score faulty')
        end
        
        if ~isnan(cum_score{iCell}(iTrg))
            last_score      = cum_score{iCell}(iTrg);
        end
    end
end

out                         = in;
out.cum_score               = in.trg_score;
out.trg_score               = new_score';
out                         = out(:,[1:17 21 18:20]);
end

function check_output(out)
%%% Check NaNs
if sum(sum(sum(isnan(out)))) > 0
    error('NaN > 0')
end

%%% Check direction
cc = 0;
for iRow = [6 8]
    cc = cc+1;
    dflag(cc) = sum(sum(out(iRow,:,:) <0)) ~= 0 || sum(sum(out(iRow,:,:) > 360)) ~= 0;
end

if sum(dflag) > 0
    error('Direction exceeds range')
end

%%% Check scores
for iRow = 11:14
    sflag(iRow-8) = sum(sum(out(iRow,:,:) < 0)) ~= 0 || sum(sum(out(iRow,:,:) > 1)) ~= 0;
end

if sum(sflag) > 0
    error('Faulty scores')
end
end

function [resultant, res_length] = integrateDots(dot_in,iState)

xIdx           	= logical(mod([1:size(dot_in{iState}{1},2)],2));            % Index for x/y position

for iFrme = 1:size(dot_in{iState},2)
    clear xpos_last ypos_last xpos ypos vs ve dist dt dot_idx mdf         	% Clear temporary variables
    
    % If first state of experiment or no dots found...
    if (iFrme == 1 && iState == 1) || sum(isnan(dot_in{iState}{iFrme}))
        resultant(iFrme)    = nan;
        res_length(iFrme)   = nan;
        continue
    end

    % If Frame 1 of stimulus state...
    if iFrme == 1 && iState > 1
        xpos_last           = dot_in{iState-1}{end}(xIdx);               	% Last frame of previous state: x-position
        ypos_last           = dot_in{iState-1}{end}(~xIdx);              	% Last frame of previous state: y-position
    elseif iFrme > 1 && ~sum(isnan(dot_in{iState}{iFrme-1}))
        xpos_last           = dot_in{iState}{iFrme-1}(xIdx);              	% Last frame: x-position
        ypos_last           = dot_in{iState}{iFrme-1}(~xIdx);              	% Last frame: y-position
    else
        resultant(iFrme)    = nan;
        res_length(iFrme)   = nan;
        continue
    end
    
    xpos                = dot_in{iState}{iFrme}(xIdx);                     	% This frame: x-position
    ypos                = dot_in{iState}{iFrme}(~xIdx);                   	% This frame: y-position
    
    for iDot = 1:length(ypos)
        % Dot position
        vs(iDot,:)      = [xpos_last(iDot), ypos_last(iDot)];             	% last frame
        ve(iDot,:)    	= [xpos(iDot), ypos(iDot)];                       	% this frame
        
        % Dot speed
        dt(iDot,:)      = ve(iDot,:) - vs(iDot,:);                         	% delta
        dist(iDot)      = pdist([vs(iDot,:);ve(iDot,:)],'euclidean');     	% distance/vector length between points
    end
    
    ofs                 = .001;                                             % Arbitrary offset
    expected_speed      = (8/1000)*(1000/120);                              % Dot speed: 8dva/s; Frame rate: 120Hz = 8.333ms/frame
    
    if median(dist) < expected_speed-ofs || median(dist) > expected_speed+ofs
%         iFrme
%         median(dist)
        resultant(iFrme) = nan;
        res_length(iFrme) = nan;
        continue
    end
    
    % Index: Exclude dots that reappear elsewhere after lifetime expired
    dot_idx            	= dist >= expected_speed-ofs & dist <= expected_speed+ofs;

    %                 % Calculate direction (theta) and step size (rho => dot speed)
    %                 [theta,rho]         = cart2pol(dt(dot_idx,1),dt(dot_idx,2));
    %                 % Zero on top, clockwise rotation
    %                 theta               = -rad2deg(theta - (pi/2));
    
    % Calculate coherence by adding up vector
    mdf                	= mean(dt(dot_idx,:));                             	% Mean x/y of all dots that didn't jump
    resultant(iFrme) 	= mod(atan2d(mdf(1),mdf(2)),360);                   % Resultant vector direction [deg]
    res_length(iFrme)	= pdist([[0 0];mdf],'euclidean') ./ median(dist);   % Resultant vector length
    end
end
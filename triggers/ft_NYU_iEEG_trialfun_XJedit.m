%% Function to read events & create file to epoch NYU iEEG data
% Function to epoch NYU iEEG data. 
% 
% fcfg.dataset = path to continuous data structure you will be epoching
% fcfg.dsfact = rate that the samples must be divided by (if applicable) - EJK need to add
% fcfg.minduration = smallest number of time between evenst (in seconds) 
% fcfg.timelim = 2 vectors of the time the events should be between (in seconds)
% fcfg.evt = vector of events to include in the final list
% fcfg.pre = time (in s) to include before the event
% fcfg.pos = time (in s) to include after the event
% fcfg.channel = {'Pulse*'};
% Created by Erik Kaestner (4-4-14) ekaestne@ucsd.edu

% Xiaojing did minor editing to make it more readable to her
% XJ modified 07/23/14, can now also be used to check if the rest-period data have been correctly aligned


function triggers = ft_NYU_iEEG_trialfun_XJedit(params, match_dur)

%% Get trigger data

cfg = [];
cfg.continuous = 'yes';
cfg.channel  = params.trigchan;

clear dat
for r = 1:params.filenum
    cfg.dataset = params.filenames{r};
    dat{r} = ft_preprocessing(cfg);
end

% For checking if the rest-period data have been correctly aligned
if nargin==2
    for clin = 1:2
        fileid = find(params.clinsys==clin);
        dat{fileid}.trial{1} = dat{fileid}.trial{1}(:, match_dur(clin,1):match_dur(clin,2));
        dat{fileid}.time{1}  = dat{fileid}.time{1}(1:length(match_dur(clin,1):match_dur(clin,2)));
    end
end
% For checking if the rest-period data have been correctly aligned


%% Set up to find the timings

Fs = dat{1}.hdr.Fs;

thresh = max(max(dat{1}.trial{1}))/5;  % about 3 times as high as the threshold used in kb_trig_mult.m
nsamp = numel(dat{1}.time{1});
ntrigchan = numel(params.trigchan);
delta   = round(0.01*Fs/2); % make sure detections for the same pulse are grouped before conversion

%%
clear triggers
for r = 1:params.filenum
    
    %% Find the timings
    
    [crsind{1:ntrigchan}] = deal([]);
    
    for k = 1:ntrigchan
        
        tmpdat  = dat{r}.trial{1}(k,:);
        ind     = crossing(tmpdat,[],thresh);
        % only keep left edge of pulse
        sel     = (tmpdat(ind+1)-tmpdat(ind-1)) > 0;
        ind     = ind(sel);
        
        crsind{k} = ind;
        clear ind sel tmpdat
    end
    
    bincode = zeros(ntrigchan,nsamp);
    for k = 1:ntrigchan 
        bincode(k,crsind{k}) = 1; 
    end
    detind  = any(bincode,1); % whether there's a detection in any trigchan    
    ind     = find(detind); % get indices at which a detection occurs in any trigchan
    
    % find detections with a second detection > delta steps before it
    tmp     = cellfun(@(x)any((ind>x-delta & (ind<x))),num2cell(ind));
    R_ind   = find(tmp);
    
    % find detections with a second detection < delta steps after it
    tmp     = cellfun(@(x)any((ind<x+delta & (ind>x))),num2cell(ind));
    L_ind   = find(tmp);
    
    % remove detections that are between two other detections in the same pulse
    Rsel    = ~ismember(R_ind,L_ind); % ismember(A,B) returns an array of the same size as A containing true where the elements of A are in B
    Lsel    = ~ismember(L_ind,R_ind);
    R_ind   = R_ind(Rsel); clear Rsel
    L_ind   = L_ind(Lsel); clear Lsel
    
    %% for each pair (L,R), set [ind(L):ind(R)] to 1 in each bincode row w/ a detection
    
    for k = 1:ntrigchan
        sel = cellfun(@(x,y)any(bincode(k,ind(x):ind(y))),num2cell(L_ind),num2cell(R_ind));
        sel = cellfun(@(x,y)ind(x):ind(y),num2cell(L_ind(sel)),num2cell(R_ind(sel)),'uniformoutput',false);
        sel = unique([sel{:}]);
        bincode(k,sel) = 1;
    end
    
    detind  = any(bincode,1); % whether there's a detection in any trigchan    
    samp    = find(detind(1:end-1)==0 & detind(2:end)==1) + 1; % only keep the first pulse in contiguous detections
        
    %% Find code for event
    
    bincode = flipud(bincode(:,samp));
    bincode = mat2cell(bincode,ntrigchan,ones(1,length(samp)));
    
    evt    = cellfun(@(x)polyval(x,2),bincode);  % convert binary to decimal
    
    triggers{r} = [samp' evt'];
    
end









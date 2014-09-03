function triggers = kb_trig_mult_XJ_edit(params, match_dur)

%% Find triggers in multiple files
% Xiaojing did minor editing to make it more readable to her
% XJ modified 07/14/14, can now detect 255 triggers - in principle any trigger codes
% XJ modified 07/23/14, can now also be used to check if the rest-period data have been correctly aligned

cfg = [];
cfg.continuous = 'yes';
cfg.channel  = params.trigchan;

clear trig_data
for r = 1:params.filenum   
    cfg.dataset = params.filenames{r};    
    trig_data{r} = ft_preprocessing(cfg);
end

% For checking if the rest-period data have been correctly aligned
if nargin==2
    for clin = 1:2
        fileid = find(params.clinsys==clin);
        trig_data{fileid}.trial{1} = trig_data{fileid}.trial{1}(:, match_dur(clin,1):match_dur(clin,2));
        trig_data{fileid}.time{1}  = trig_data{fileid}.time{1}(1:length(match_dur(clin,1):match_dur(clin,2)));
    end
end
% For checking if the rest-period data have been correctly aligned

clear triggers
for r = 1:params.filenum
    
    my_trl_dc = [];
    
    for itrig = 1:8
        
        dc = find(trig_data{r}.trial{1}(itrig,:) > params.trig_thresh); % default = 190000, samples above threshold; will later refer back to this matrix
        dc_sub = diff(dc); % subtract sample with preceding sample (e.g., 1 1 1 1 368 1 1 1 1 600 1 1 1 1)
        
        if length(dc_sub) >2
            
            dc_sub_shift = [dc_sub, 1]; % add one to end
            dc_sub_first = [600,dc_sub]; % add first value of 600 to dc_sub, so it is recognized as trigger (e.g., 360 1 1 1 1 368 1 1 1 1 600 1 1 1 1)
            dc_event = find(dc_sub_first > 10 & dc_sub_shift == 1); % find samples with big gaps to the pervious trigger, with also the immediate next sample > params.trig_thresh
            trigz_dc = dc(dc_event)';
            %dc4_trl = [trigz_dc4 - 250, trigz_dc4 + 250, repmat(-250,length(trigz_dc4),1)];

            my_trl_dc = [my_trl_dc; [trigz_dc, repmat(itrig,length(trigz_dc),1)]];

        end
    end
    
    triggers{r} = sortrows(my_trl_dc,1);
    triggers{r}(~all(triggers{r},2),:) = [];
    
end


for r = 1:params.filenum
    
    triggers{r}(1,4) = [0]; %creates a 4th col
    triggers{r}(:,3) = pow2(triggers{r}(:,2)-1);
    
    ieve = 1;
    seed = 0;  % Counting the number of rows belonging to the same trigger
    while ieve < size(triggers{r},1);
        if abs(triggers{r}(ieve+1,1)-triggers{r}(ieve,1)) < 10;
            if seed == 0
                triggers{r}(ieve,4) = triggers{r}(ieve,3) + triggers{r}(ieve+1,3);
            elseif seed > 0
                triggers{r}(ieve,4) = triggers{r}(ieve,4) + triggers{r}(ieve+1,3);
            end
            seed = seed + 1;
            triggers{r}(ieve+1,:) = [];
        else
            if seed > 0
                seed = 0;
            end
            if triggers{r}(ieve,4) == 0;
                triggers{r}(ieve,4) = triggers{r}(ieve,3);
            end
            ieve = ieve + 1;
        end
        
    end
    
    if  triggers{r}(end,4) == 0;
        triggers{r}(end,4) = triggers{r}(end,3);
    end
    
    triggers{r}(:,2:3) = [];
    
end




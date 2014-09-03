
function varargout = good_EEG_channels_fn(params, cont_data)


%% Figure out which files are clin1/2

clear clinsys
for r = 1:params.filenum
    if(strfind(lower(params.filenames{r}),'clin1'))
        clinsys(r) = 1;
    elseif(strfind(lower(params.filenames{r}),'clin2'))
        clinsys(r) = 2;
    end 
end
varargout{1} = clinsys;


if nargin==2
    
    clear EEG_channels_all
    for clin = 1:2
        
        %% Find out if all the data files have the same channels with the same order
        
        fileid = find(clinsys==clin);
        
        Ndata = length(fileid);
        Nchan = zeros(1,Ndata);
        label = {};
        for i=1:Ndata
            Nchan(i)  = length(cont_data{fileid(i)}.label);
            label = cat(1, label(:), cont_data{fileid(i)}.label(:));
        end
        
        alllabel = unique(label, 'first');
        order    = zeros(length(alllabel),Ndata);
        for j=1:Ndata
            tmplabel = cont_data{fileid(j)}.label;
            [ix,iy]  = match_str(alllabel, tmplabel);  % This is a fieldtrip function
            order(ix,j) = iy;
        end
        
        cmp_label = ~all(all(order-repmat(order(:,1),[1 Ndata])==0));
        if(cmp_label)
            disp('Not all data files have the same channels and/or ordering!');
            pause;
        end
        
        %% Select all EEG channels
        
        all_channels = cont_data{fileid(1)}.label;
        EEG_channels = {};
        for ch = 1:length(all_channels)
            if(strcmp(all_channels{ch}(1:3), 'EEG'))
                EEG_channels = [EEG_channels; all_channels(ch)];
            end
        end
        
        %% Remove bad channels identified using ft_databrowser
        
        bad_chs = params.bad_chs{clin};
        
        bad_ch_idx = [];
        for ch=1:length(bad_chs)
            bad_ch_name = ['EEG ' bad_chs{ch} '-REF'];
            for ach=1:length(EEG_channels)
                if(strcmp(EEG_channels{ach},bad_ch_name))
                    bad_ch_idx(ch) = ach;
                end
            end
        end
        EEG_channels(bad_ch_idx)=[];
        
        EEG_channels_all{clin} = EEG_channels;
        
    end
    
    varargout{2} = EEG_channels_all;
end








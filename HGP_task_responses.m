% Integrating Kristen's & UCSD scripts


%% Set parameters

%Initialize fieldtrip, if this errors ensure that fieldtrip is in path
ft_defaults;

%Calls on meta_file to re-initialize params, dataset_idx contains
%data-specific values
dataset_idx = 4;
params = meta_file_july28(dataset_idx); % The function clears params

%Frequency bands of interest for analysis
params.linenoise_freq = [60 120 180]; % Hz
params.hgp_bpfreq = [70 190]; % Hz

if(params.RunOrRest)
    params.blc_window = [-0.4 0]; % in sec, Baseline correction window for epoched data
    params.rejvis_lat = [-0.5 0.8]; % in sec, ft_rejectvisual latency
end

%% Make analysis (plot) directory, load raw data

analysis_dir = [params.main_dir filesep params.experiment filesep params.subject filesep 'PD']; 
if exist(analysis_dir,'dir')==0
    mkdir(analysis_dir);
end
params.analysis_dir = analysis_dir;

cfg = [];
cfg.continuous = 'yes';
cont_data = {};
for i = 1:params.filenum
    cfg.dataset = params.filenames{i};
    cont_data{i} = ft_preprocessing(cfg);
end

%% Plot raw data, identify very bad channels

file_num = 3;
cfg = [];
cfg.continuous = 'yes';
cfg.viewmode = 'vertical';

ft_databrowser(cfg,cont_data{file_num});

%% Detect trigger times and codes

triggers_method1 = kb_trig_mult_XJ_edit(params); % Adapted from Kristen's script
triggers_method2 = ft_NYU_iEEG_trialfun_XJedit(params); % Adapted from UCSD script 

% compare trigger times & codes derived by the two methods
for r = 1:params.filenum
    
    T1 = triggers_method1{r};
    T2 = triggers_method2{r};
    
    if size(T1,1)~=size(T2,1) || sum(abs(T1(:,2)-T2(:,2)))~=0
        error('Major error with trigger detection!');
    else
        mean(abs(T1(:,1)-T2(:,1)))
    end
end

if strcmpi(params.experiment,'PM')
    [triggers_method3] = PM_addresponsetrigs(params.logfile,triggers_method2,params.events);
    params.my_trl = triggers_method3;

else
    params.my_trl = triggers_method2; % Use UCSD triggers since they are earlier by 1 sample...
end
%% Define trials

if(params.RunOrRest)
    
    Fs = cont_data{1}.hdr.Fs;
    pre_len = 1;  % before trigger onset, in sec
    post_len = 2;  % after trigger onset, in sec
    params.pre_samp = round(pre_len*Fs);
    params.post_samp = round(post_len*Fs);
    
    clear params.trials
    for r = 1:params.filenum
        params.trials{r} = [params.my_trl{r}(:,1) - params.pre_samp,...
            params.my_trl{r}(:,1) + params.post_samp,...
            repmat(-params.pre_samp, size(params.my_trl{r},1), 1),...
            params.my_trl{r}(:,2)];
    end
    
end

%% Select good EEG channels

[params.clinsys, params.EEG_goodchan_perclin]...
                    = good_EEG_channels_fn(params,cont_data); % 1. Find out if all the data files have the same channels with the same order
                                                              % 2. Select all EEG channels
                                                              % 3. Remove bad channels identified using ft_databrowser

%% Compute epoched HGP data for tasks, or continuous HGP data for rest periods
%hgp_data_clin is trial data separated between both clinical systems
%hgp_trials is trial data with the clinical systems combined
hgp_data_clin = cell(1,2);

for clin = 1:2
    
    fileid = find(params.clinsys==clin);
    
    cont_data = {};
    cont_data_hgp = {};
    epoch_data_hgp = {};
    
    for fl = 1:length(fileid)
        
        % Load raw data from selected EEG channels
        cfg = [];
        cfg.continuous = 'yes';
        cfg.dataset = params.filenames{fileid(fl)};
        cfg.channel = params.EEG_goodchan_perclin{clin};
        cont_data{fl} = ft_preprocessing(cfg);
        
        % Get rid of line noise
        cfg = [];
        cfg.dftfilter = 'yes';
        cfg.dftfreq = params.linenoise_freq;
        cont_data{fl} = ft_preprocessing(cfg, cont_data{fl});
        
        % Apply average referencing
        cfg = [];
        cfg.reref = 'yes';
        cfg.refchannel = params.EEG_goodchan_perclin{clin};   % let data = cont_data{fl}.trial{1};
        cont_data{fl} = ft_preprocessing(cfg, cont_data{fl}); % this is equivalent to --> data - repmat(mean(data,1), size(data,1), 1)
        
        % Apply HGP filter
        cfg = [];
        cfg.bpfilter = 'yes';
        cfg.bpfreq = params.hgp_bpfreq;
        cont_data_hgp{fl} = ft_preprocessing(cfg, cont_data{fl});
        
        % Apply Hilbert transformation
        cfg = [];
        cfg.hilbert = 'abs';
        cont_data_hgp{fl} = ft_preprocessing(cfg, cont_data_hgp{fl});
        
        % Epoch the HGP data
        if(params.RunOrRest)
            cfg = [];
            cfg.trl = params.trials{fileid(fl)};
            epoch_data_hgp{fl} = ft_redefinetrial(cfg, cont_data_hgp{fl});
        end
    end
     
    % Baseline correction
        cfg = [];
        cfg.demean = 'yes';
        cfg.baselinewindow = params.blc_window;
        
    if(params.RunOrRest)
        %Fix this
        % Append different runs of the same clinical system
%         if length(epoch_data_hgp{clin}) == 1
            cfg = [];
            hgp_data_clin{clin} = epoch_data_hgp{1};
%         else
%             cfg = [];
%             hgp_data_clin{clin} = ft_appenddata(cfg,epoch_data_hgp{1:end});
%         end
        %?
        %hgp_data_clin{clin} = ft_preprocessing(cfg, epoch_data_hgp{clin});  % All samples in each single trial subtracts a constant
    else
        if length(cont_data_hgp) == 1
            hgp_data_clin{clin} = cont_data_hgp{1};
        end
    end
        
end

% Append the two clinical systems
 cfg = [];
 hgp_trials = ft_appenddata(cfg, hgp_data_clin{1:end});  % {1:end} is necessary!!


%% Remove noisy trials/channels

if(params.RunOrRest) % rest-period recordings do not have trials - skip ft_rejectvisual
    
    cfg = [];
    cfg.method = 'summary';
    cfg.latency = params.rejvis_lat;
    
    hgp_trials_final = ft_rejectvisual(cfg, hgp_trials);
else
    hgp_trials_final = hgp_trials;
end

params.EEG_goodchan_all = hgp_trials_final.label;


%% Separate trials by event type, calculate mean and SEM across trials of the same event type
%Creates hgp_trials_avgsem to use for plotting

dat_in_mat = cat(3, hgp_trials_final.trial{:});

clear hgp_trials_avgsem

time_in_mat = squeeze(cat(3, hgp_trials_final.time{:}))'; % This won't work if not all trials have the same number of samples

if(all(all((time_in_mat - repmat(hgp_trials_final.time{1}, size(time_in_mat,1), 1))==0)))
    hgp_trials_avgsem.timeavg = hgp_trials_final.time{1};
else
    error('Error! Not all trials have the same time points!');
end

[ieve_idx, hgp_trials_avgsem.avgdat] = deal(cell(1, length(params.events)));

for ieve = 1:length(params.events)
    
    ieve_idx{ieve} = find(hgp_trials_final.trialinfo == params.events(ieve));
    
    hgp_trials_avgsem.avgdat{ieve}.avg = mean(dat_in_mat(:,:,ieve_idx{ieve}), 3);
    hgp_trials_avgsem.avgdat{ieve}.sem = std(dat_in_mat(:,:,ieve_idx{ieve}), [], 3)/sqrt(length(ieve_idx{ieve}));
                     %%%%%%%%%%% Note %%%%%%%%%%%%
                     % ft_timelockanalysis uses the 'Na?ve algorithm' for
                     % calculating variance, see http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
                     % which is NOT recommended
                     % To convert the var calculated by ft_timelockanalysis
                     % into SEM, do
                     % sqrt(avg_data{ieve}.var(ichan,:)) ./ sqrt(length(avg_data{ieve}.trialinfo))
                     %%%%%%%%%%% Note %%%%%%%%%%%%
    hgp_trials_avgsem.avgdat{ieve}.trialinfo = hgp_trials_final.trialinfo(ieve_idx{ieve});
    
end


%% Plot HGP

% params.plotchan = 'all';
params.prestim_plot = -0.3; % in sec, for plotting only
params.poststim_plot = 1; % in sec, for plotting only
params.plotevt = params.events;

visible = 1; % =1 necessary for printing figures with the correct dimensions

plot_avgHGP_allchan_fn_PD(hgp_trials_avgsem.avgdat{1:3}, params, visible);






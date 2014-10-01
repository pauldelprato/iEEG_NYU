   % Integrating Kristen's & UCSD scripts


%% Set parameters

%Initialize fieldtrip, if this errors ensure that fieldtrip is in path
ft_defaults;

%Calls on meta_file to re-initialize params, dataset_idx contains
%data-specific values
dataset_idx = 4;
params = meta_file_v4(dataset_idx); % The function clears params

%Frequency bands of interest for analysis
params.linenoise_freq = [60 120 180]; % Hz
params.hgp_bpfreq = [70 190]; % Hz

if(params.RunOrRest)
    params.blc_window = [-.2 0]; % in sec, Baseline correction window for epoched data
    params.rejvis_lat = [0 10]; % in sec, ft_rejectvisual latency
end

params.Gaussian_smooth = 1;
params.gausssd = 10;
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
    pre_len = .2;  % before trigger onset, in sec
    post_len = 10;  % after trigger onset, in sec
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

if params.Gaussian_smooth
    Fs = cont_data{1}.hdr.Fs;
    gausssd = params.gausssd/1000*Fs;
    
    gauss_x = -20:20;
    gauss_y = normpdf(gauss_x,0,gausssd);
    gauss = gauss_y/sum(gauss_y);
end

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
        
        %Apply Gaussian smoothing -- optional
        if params.Gaussian_smooth
            cont_data_hgp{fl}.trial{1}=convn(cont_data_hgp{fl}.trial{1},gauss,'same');
        end
        
        % Epoch the HGP data
        if(params.RunOrRest)
            cfg = [];
            cfg.trl = params.trials{fileid(fl)};
            epoch_data_hgp{fl} = ft_redefinetrial(cfg, cont_data_hgp{fl});
        end
    end
     
   
        
    if(params.RunOrRest)
        %Fix this
        % Append different runs of the same clinical system
        if length(epoch_data_hgp) == 1
%             cfg = [];
%             hgp_data_clin{clin} = epoch_data_hgp;
            
           % Baseline correction
        cfg = [];
        cfg.demean = 'yes';
        cfg.baselinewindow = params.blc_window;
        hgp_data_clin{clin} = ft_preprocessing(cfg,epoch_data_hgp{1});
        else
            hgp_data_clin{clin} = ft_appenddata(cfg,epoch_data_hgp{1:end});
            cfg = [];
        cfg.demean = 'yes';
        cfg.baselinewindow = params.blc_window;
        hgp_data_clin{clin} = ft_preprocessing(cfg,epoch_data_hgp{clin});
        end
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


%%
% Determine pre-stimulus relevant samples for PM

prelen_samp = zeros(1,length(RTs_outrem));

for i = 1 : length(RTs_outrem)
    prelen_samp(i) = fix(RTs_outrem(i) * (512000/1000)); %Converts from tens of ms to samples rounded down
end

% for i = 1: length(prelen_samp)
%     hgp_trials_final_rt_fixed = hgp_trials_final.trial{i}(prelen_samp

%%
params.combinedevents = {[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]};
hgp_trials_comb = separateEventdata(hgp_trials_final,params.combinedevents);
%% 
%Plot HGP
test_electrode = 'GB_54';%idx 102
test_electrode2 = 'GA_64';%idx 48

 params.plotchan = params.buttonpress_electrodes;
params.prestim_plot = -.2; % in sec, for plotting only
params.poststim_plot = 9; % in sec, for plotting only
params.plotevt = params.combinedevents;
plot_avgHGP_wide(hgp_trials_comb, params,params.plotchan(1));
%Aud stim
% blc_begin = sec_to_sample(0);
% blc_end = sec_to_sample(.2);
% time_end = sec_to_sample(.78+.2);
% [hgp_trials_comb2.avgdat{1}.avg(:,blc_begin:time_end)] = ft_preproc_polyremoval(hgp_trials_comb.avgdat{1}.avg(:,blc_begin:time_end), 0, blc_begin, blc_end)
% %  params.plotchan = 'all';
% % params.prestim_plot = -.2; % in sec, for plotting only
% % params.poststim_plot = .75; % in sec, for plotting only
% params.plotevt = params.combinedevents;
% plot_avgHGP(hgp_trials_comb2, params);
%plot_avgHGP(hgp_trials_avgsem, params,params.auditory_electrodes);
params.auditory_electrodes = {'GA_07','GA_08','GA_14','GA_15','GA_20','GA_28','GB_39',...
    'GB_40','GB_45','GB_46','GB_47','GB_48','GB_53','GB_54','GB_55','GB_56','GB_61',...
    'GB_62','GB_63','GB_64','DPI_01','DPI_02','DPI_03'};

%Vis stim
%relative to data input
% blc_begin = sec_to_sample(.78+.2);
% blc_end = sec_to_sample(.8+.2);
% time_end = sec_to_sample(1.5+.2);
% blc_rel_begin = sec_to_sample(0);
% blc_rel_end = sec_to_sample(.2);
% [hgp_trials_comb2.avgdat{1}.avg(:,blc_begin:time_end)] = ft_preproc_polyremoval(hgp_trials_comb.avgdat{1}.avg(:,blc_begin:time_end), 0, blc_rel_begin, blc_rel_end)
%  params.plotchan = 'all';
% params.prestim_plot = .78; % in sec, for plotting only
% params.poststim_plot = 1.5; % in sec, for plotting only
% params.plotevt = params.combinedevents;
% 
% visible = 1; % =1 necessary for printing figures with the correct dimensions
% plot_avgHGP(hgp_trials_comb2, params);
params.visual_electrodes = {'GA_14','MT_01','PT_04','O_05','O_06'};
%Responses
% blc_begin = sec_to_sample(1.5+.2);
% blc_end = sec_to_sample(1.8+.2);
% time_end = sec_to_sample(9+.2);
% blc_rel_begin = sec_to_sample(0);
% blc_rel_end = sec_to_sample(.3);
% [hgp_trials_comb2.avgdat{1}.avg(:,blc_begin:time_end)] = ft_preproc_polyremoval(hgp_trials_comb.avgdat{1}.avg(:,blc_begin:time_end), 0, blc_rel_begin, blc_rel_end)
% %  params.plotchan = 'all';
% % params.prestim_plot = 1.5; % in sec, for plotting only
% % params.poststim_plot = 9; % in sec, for plotting only
% params.plotevt = params.combinedevents;
% 
% visible = 1; % =1 necessary for printing figures with the correct dimensions
% plot_avgHGP(hgp_trials_comb2, params);
params.buttonpress_electrodes = {'GA_49','GA_63','GA_64','GB_04','GB_05','GB_54','GB_55'};

%% Compute epoched HGP data for tasks, or continuous HGP data for rest periods
% !!! NO BLC !!!


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
     
   
        
    if(params.RunOrRest)
        %Fix this
        % Append different runs of the same clinical system
        if length(epoch_data_hgp) == 1
%             cfg = [];
%             hgp_data_clin{clin} = epoch_data_hgp;
            
           % Baseline correction
        cfg = [];
%         cfg.demean = 'yes';
%         cfg.baselinewindow = params.blc_window;
        hgp_data_clin_no_blc{clin} = ft_preprocessing(cfg,epoch_data_hgp{1});
        else
            hgp_data_clin{clin} = ft_appenddata(cfg,epoch_data_hgp{1:end});
            cfg = [];
%         cfg.demean = 'yes';
%         cfg.baselinewindow = params.blc_window;
        hgp_data_clin_no_blc{clin} = ft_preprocessing(cfg,epoch_data_hgp{clin});
        end
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
if params.RunOrRest
% Append the two clinical systems
 cfg = [];
 hgp_trials = ft_appenddata(cfg, hgp_data_clin_no_blc{1:end});  % {1:end} is necessary!!
else
cfg = [];
 hgp_trials = ft_appenddata(cfg, hgp_data_clin{1:end});  % {1:end} is necessary!!
end

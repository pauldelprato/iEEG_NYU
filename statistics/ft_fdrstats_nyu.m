%% Function to calculate either a ttest or anova and apply FDR corrections
%
% ft_data
%
% Outputs:
% data: fieldtrip stat data structure % - EJK need to find out what a ft_stats structure looks like and modify
%                                     % script to return that
%
% Parameters: default : options :  description
% cfg.events     : []  : [] : event code to compare, note One-Way ANOVA will only be performed on all specified event codes
% cfg.alpha'     : .05 : [] :  alpha level to be specified used for the  FDR correction, the default value for alpha is set at 5% significance level
%                             ANOVA is performed based on 5% significance as a matlab built in defult for anova1() function
%                             Assumptions: All sample populations are normally distributed
%                                          All sample populations have equal variance
%                                          All observations are mutually independent
% cfg.fdr_type   :'row': {'row','col','mtx'} : specification on conducting FDR computations:
%                       Assumptions: This FDR method assumes indepence within passed in p-values (and their corresponding F-statistic)
%                       'mtx' would do computations on the entire matrix as a whole, (therefore the FDR results for this option would be the most conservative)
%                       'col' would do computations across all the channels at each time point
%                       'row' would do computations across all time points in each channel
% cfg.basetime : If specified, will calculate baseline vs chosen events,
%                       baseline defined as [start end] in seconds
%
%
% created by Nadira Babaeva on 2013/09/06
% edited on 2013/10/24
% edited on 2014/4/25 by Erik Kaestner (ekaestne@ucsd.edu) to change to
%       accept fieldtrip inputs & to calculate baseline stats


function [stat_data] = ft_fdrstats_nyu(cfg,ft_data)
%% handle inputs
% sanity check
if isempty(cfg.events) || (numel(cfg.events) < 2 && ~isfield(cfg,'basetime'))
    error('data must have at least two events.');
elseif numel(cfg.events) > 1 && isfield(cfg,'basetime')
    error('can only run stats on one baseline at a time');
end

if ~isfield(cfg,'fdr_type')
    cfg.fdr_type = 'row';
end
if ~isfield(cfg,'alpha')
    cfg.alpha = .05;
end


%% One-Way Anova (Unbalanced design)
ft_dat = cat(3,ft_data.trial{:});
if isfield(cfg,'basetime')
    bcfg = [];
    bcfg.basetime = [cfg.basetime(1) cfg.basetime(2)];
    base_dat = ft_base_create(bcfg,ft_data);
    base_dat = cat(3,base_dat.trial{:});
end

% finding the max number of trials throughout all conditions:
sizevect = zeros(1,numel(unique(cfg.events)));
ev_ind   = cell(1,numel(unique(cfg.events)));
for ievent=1:numel(sizevect)
    ev_ind{ievent} = find(ft_data.trialinfo==cfg.events(ievent));
    sizevect(ievent) = numel(ev_ind{ievent});
end
maxtrials = max(sizevect);
mtrx = zeros([maxtrials numel(cfg.events)]);

% Creating an empty matrix to store p-values
Pmtx   = zeros([numel(ft_data.label) length(ft_data.time{1})]);
Fstats = cell([numel(ft_data.label) length(ft_data.time{1})]);
clear Fstats
% Performing One-Way ANOVA on all the trials in each contion(as gropus)
% for each channel at each point in time:
for ichan=1:numel(ft_data.label)
    for iti=1:length(ft_data.time{1})
        for ievent=1:numel(cfg.events)
            depthvect = squeeze(ft_dat(ichan,iti,cell2mat(ev_ind(ievent))));
            depthvect = cat(1,depthvect,NaN(maxtrials-sizevect(ievent),1));
            mtrx(:,ievent) = depthvect;
            if exist('base_dat','var')
                depthvect = squeeze(base_dat(ichan,iti,cell2mat(ev_ind(ievent))));
                depthvect = cat(1,depthvect,NaN(maxtrials-sizevect(ievent),1));
                mtrx(:,ievent+1) = depthvect;
            end
        end
        [Pmtx(ichan,iti),~,Fstats(ichan,iti)]=anova1(mtrx, [], 'off');
    end
end

%mmil_logstr('Finished running One-Way Anova')
disp('Finished running One-Way Anova')
%% Perform FDR

% Computing the False Discovery Rate (FDR) matrix here:
switch cfg.fdr_type
    case 'mtx'
        [~,fdrmask] = fdr_cor([Pmtx(:)]',cfg.alpha);
        fdrmask = vec2mat(fdrmask,size(Pmtx,1))';
    case 'col'
        fdrmask=zeros(size(Pmtx))';
        for itim = 1:size(Pmtx,2)
            [~, fdrmask(itim,:)]=fdr_cor([Pmtx(:,itim)]',cfg.alpha);
        end
        fdrmask = fdrmask';
    case 'row'
        fdrmask=zeros(size(Pmtx));
        for ichan = 1:size(Pmtx,1)
            [~, fdrmask(ichan,:)]=fdr_cor(Pmtx(ichan,:),cfg.alpha);
        end
end

disp('Finished running False Discovery Rate corrections')
%mmil_logstr('Finished running False Discovery Rate corrections')

% Setting invalid P-values to 1 (applying the FDR corrections)
Pmtx(~fdrmask) = 1;
disp('All False Discovery Rate corrections have been applied')
%mmil_logstr('All False Discovery Rate corrections have been applied')

% General structure of stat_data:
stat_data.num_sensors = double(size(Pmtx,1));
stat_data.sfreq = double(size(Pmtx,2));
stat_data.sensor_info = [ft_data.label];
stat_data.stats = struct;
stat_data.cfg = struct;

%stat_data = ft_addscript_nyu(stat_data);

% Subtructures of stats field in stat_data for event code -1:
stat_data.stats(1).mask = zeros(size(Pmtx));
stat_data.stats(1).prob = Pmtx;
stat_data.stats(1).stat = [vec2mat([Fstats.s], size(Pmtx,1))]';
stat_data.stats(1).dimord = 'chan_time';
stat_data.stats(1).posclusters = [];
%  Subfields of posclusters:
stat_data.stats(1).posclusters.prob = [];
stat_data.stats(1).posclusters.clusterstat = [];
stat_data.stats(1).negclusters = [];
%  Subfields of negclusters:
stat_data.stats(1).negclusters.prob = [];
stat_data.stats(1).negclusters.clusterstat = [];

stat_data.stats(1).negclusterslabelmat = struct;
stat_data.stats(1).posclusterslabelmat = struct;
stat_data.stats(1).time = [ft_data.time];
stat_data.stats(1).label = {ft_data.label};
stat_data.stats(1).cfg = [];
stat_data.stats(1).event_code = -1;

for ievent = 1:numel(cfg.events)
    % Subtructures of stats field in stat_data for remaining event codes:
    stat_data.stats(ievent+1) = stat_data.stats(1);
    stat_data.stats(ievent+1).event_code = cfg.events(ievent);
end

% Subtructures of params field in stat_data:
stat_data.cfg.events = [cfg.events];

% FDR function code:
    function [fdrthresh,fdrmask] = fdr_cor(pval_vec,alpha)
        % Compute FDR and return the fdrthresh and a mask of significance
        %   created by Brian Quinn
        pval_sorted = sort(pval_vec);
        fdrpvec=((1:length(pval_vec))./length(pval_vec))*alpha;
        fdrtest = pval_sorted-fdrpvec;
        fdrp = find(fdrtest < 0,1,'last');
        fdrmask = zeros(1,length(pval_vec));
        if ~isempty(fdrp),
            fdrthresh = pval_sorted(fdrp);
            fdrmask(pval_vec < fdrthresh)=1;
        else
            fdrthresh=0;
        end
    end

end
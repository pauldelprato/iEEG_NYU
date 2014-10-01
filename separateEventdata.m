function [hgp_trials_avgsem] = separateEventdata(hgp_trials_final,events)

dat_in_mat = cat(3, hgp_trials_final.trial{:});

clear hgp_trials_avgsem

time_in_mat = squeeze(cat(3, hgp_trials_final.time{:}))'; % This won't work if not all trials have the same number of samples

if(all(all((time_in_mat - repmat(hgp_trials_final.time{1}, size(time_in_mat,1), 1))==0)))
    hgp_trials_avgsem.timeavg = hgp_trials_final.time{1};
else
    error('Error! Not all trials have the same time points!');
end

[ieve_idx, hgp_trials_avgsem.avgdat] = deal(cell(1, length(events)));

for ieve = 1:length(events)
    ieve_idx{ieve} = [];
    
    if iscell(events)
        for j = 1 : length(events{ieve})
        ieve_idx{ieve} = [ieve_idx{ieve} ; find(hgp_trials_final.trialinfo == events{ieve}(j))];
        end
    else
    ieve_idx{ieve} = find(hgp_trials_final.trialinfo == events(ieve));
    end
    
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
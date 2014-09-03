function [final_triggers] = PM_addresponsetrigs(filename,trigger_list,event_codes)

if iscell(trigger_list)
    num_files = length(trigger_list);
else
    num_files = 1;
    trig_input = trigger_list;
    trigger_list = {};
    trigger_list{1} = trig_input;
    
end

log_output = mmil_readtext(filename,'\t');

%remove header
log_output(1:5,:) = [];


%find all log times and trigger codes
log_triggers = [];
log_triggers_no64 = [];
for i = 1 : length(log_output)
    %trigger code is a number and is not 0 (0=visual word pres.)
    if and(isnumeric(cell2mat(log_output(i,4))),cell2mat(log_output(i,4))~=0)
        log_triggers = [log_triggers; [cell2mat(log_output(i,5)), cell2mat(log_output(i,4))]];
        %additional constraint (not the response) for additional trigger
        %log
        if cell2mat(log_output(i,4))~=64
            log_triggers_no64 = [log_triggers_no64; cell2mat(log_output(i,5)),cell2mat(log_output(i,4))];
        end
    end
end

log_sample_triggers = [ceil(log_triggers(:,1).*(512/10000)), log_triggers(:,2)];
log_sample_triggers_no64 = [log_triggers_no64(:,1).*(512/10000), log_triggers_no64(:,2)];



adjusted_triggers = {};
%read in the trigger list, determine accurate samples length of
%trial, and then add the "offsetted" samples of the response triggers
for i = 1 : num_files
    clear trigs;
    trigs = trigger_list{i};
    
    if trigs(:,2) == log_sample_triggers_no64(:,2)
        disp('All trigger codes match.');
    else
        error('All trigger codes do not match between log and trigger list.');
    end
    
    %        %hard coded trig codes
    tone = 2;
    tone_counter = 0;
    aud_word_pres = 1:15;
    aud_word_pres_counter = 0;
    
    for j = 1 : length(trigs)
        while tone_counter == 0 || aud_word_pres_counter == 0
            
            if trigs(j,2) == 2
                tone_counter=1;
                samp1_trig = trigs(j,1);
                samp1_log = log_sample_triggers(j,1);
            end
            if ismember(trigs(j,2),aud_word_pres)
                aud_word_pres_counter = 1;
                samp2_trig = trigs(j,1);
                samp2_log = log_sample_triggers(j,1);
            end
        end
    end
    
    if (samp2_trig-samp1_trig)-(samp2_log-samp1_log)<2
        disp('Log and trigger list find same trial length');
    else
        error('Different trial length found between log and trigger list.');
    end
    
    samplediff = samp1_trig-samp1_log;
    
    adjusted_triggers{i} = [log_sample_triggers(:,1)+samplediff,log_sample_triggers(:,2)];
    
    
    %remove tone trials
    
        adjusted_triggers{i}(1:7:end,:)=[];
        
        %remove non-relevant event codes
        
           final_triggers{i} = [];
        for j = 1 : length(adjusted_triggers{i})
            if ismember(adjusted_triggers{i}(j,2),event_codes)
                final_triggers{i} = [final_triggers{i}; adjusted_triggers{i}(j,:)];
            end
        end
         
       
end
end

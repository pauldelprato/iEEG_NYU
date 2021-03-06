function plot_avgHGP_allchan_fn(avg_data, params, visible)


%% Setting some variables

proc_name = 'HGP';

%Directory to save the plots
% plotdir = sprintf('%s/%s_plots',params.analysis_dir,proc_name);
% if(~exist(plotdir,'dir'))
%     mkdir(plotdir); 
% end
plotdir = params.analysis_dir;

channels = 1:length(params.EEG_goodchan_all);

line_colors = [[0 1 0];[0 0 1];[1 0 0]];

%% Events and channel number per figure

clear yevent event_name legend_names
for i = 1:length(params.plotevt)
    yevent(i) = find(params.events==params.plotevt(i));
    event_name{i} = params.event_names{yevent(i)};
    legend_names{i} = [event_name{i} ' N=' num2str(length(avg_data.avgdat{yevent(i)}.trialinfo))];
end

% Plot 25 channels per figure
chan_length = length(channels);
plot_length = floor(chan_length/25);
if chan_length > plot_length*25
    plot_length = plot_length+1;
end

xft = avg_data.timeavg;
begbin = max([nearest(xft, params.prestim_plot), 20]);  % omit first and last 20 samples because of edge effect
endbin = min([nearest(xft, params.poststim_plot),length(xft)-20]);

%% Plot

for iplot = 1:plot_length
    
    h=custom_figure(50,50,visible);
    
    if iplot == plot_length && rem(chan_length,25)~=0
        panel_num = rem(chan_length,25);
    else
        panel_num = 25;
    end
    
    for k = 1:panel_num
        
        h(k) = subplot(5,5,k);
        subplot(h(k));
        ichan = k+(iplot-1)*25;
        
        chan = params.EEG_goodchan_all{ichan};
        chan = strrep(chan,'_','-');
        
        clear yft ysem
        for ieve = 1:length(yevent)
            
            yft(ieve,:) = avg_data.avgdat{yevent(ieve)}.avg(ichan,:);
            ysem(ieve,:)= avg_data.avgdat{yevent(ieve)}.sem(ichan,:);
                     
            confplot_3andC(xft, yft(ieve,:), ysem(ieve,:), line_colors(ieve,:)); 
            hold on;
            
        end
        
        ymax_val = max(max(yft(:,begbin:endbin) + ysem(:,begbin:endbin)))*1.2; %largest y-axis value
        ymin_val = min(min(yft(:,begbin:endbin) - ysem(:,begbin:endbin)))*1.2; %smallest y-axis value
        if(ymin_val>0)
            ymin_val = 0;
        end
        hold on, plot([0 0],[ymin_val ymax_val],'color',[0.6 0.6 0.6]);
        
        axis([params.prestim_plot params.poststim_plot ymin_val ymax_val]);
        
        set(gcf, 'Color', 'w');
        title(sprintf('%s %s',proc_name,chan),'FontSize',16);
        if k == 1
            xlabel('Time (sec)');
        end
        
    end
    
    %creates legend
    lgdPosition = [0.04 0.88 0.05 0.05];
    lgdUnits = 'normalized';
    hleg = legend(legend_names);
    set(hleg, 'Box', 'off'); hold on;
    set(hleg,'Position', lgdPosition,'Units', lgdUnits);
       
    export_fig(sprintf('%s/%s_%s_%s_per_channel_%d.png',...
        plotdir,params.subject,params.experiment,proc_name,iplot),'-a1','-nocrop');
    close(gcf);
end





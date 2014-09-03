%% Function to create fakebaseline structures (which mimics the structure of the data it's the baseline of)
% To be used with ft_macro_plot & ft_NS5_plot
% 
% ft_data
% cfg.basetime = baseline time for the data (in seconds)
% 
% created by Erik Kaestner 4-25-14 (ekaestne@ucsd.edu)

function [ft_base] = ft_base_create(cfg,ft_data)

[~,bstm(1)] = min(abs(cfg.basetime(1) - ft_data.time{1}));
[~,bstm(2)] = min(abs(cfg.basetime(2) - ft_data.time{1}));

bl = cellfun(@(x) mean(x(:,bstm(1):bstm(2)),2),ft_data.trial,'UniformOutput',0);

for it = 1:numel(ft_data.trial)
    
    ft_base.trial{it} = repmat(bl{it},1,numel(ft_data.time{1}));
    
end

end
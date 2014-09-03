function h=custom_figure(width,height,visible)
%function h=custom_figure(width,height)
% Customized figure function. Opens a figure with 
% specified width and height (in centimeters)
% The properties are set such that Matlab won't resize the 
% figure while printing or exporting to graphics files (eps, 
% tiff, jpeg, ...).
%
% Get the screen size in centimeters
set(0,'units','centimeters')
scrsz=get(0,'screensize');
% Calculate the position of the figure
%position=[scrsz(3)/2-width/2 scrsz(4)/2-height/2 width height];
position=[1 2 width height];

%Suppress graphics output if desired
if visible == 0
    h=figure('Toolbar','none','visible','off');
else
    h=figure('Toolbar','none');
end

set(h,'units','centimeters')
% Place the figure
set(h,'position',position)
oposition = get(h,'OuterPosition');
set(h,'OuterPosition',[0 scrsz(4)-oposition(4) oposition(3) oposition(4)])
% Do not allow Matlab to resize the figure while printing
set(h,'paperpositionmode','auto')
% Set screen and figure units back to pixels
set(0,'units','pixel')
set(h,'units','pixel')
set(h,'Color',[1 1 1])

%
% Note: in order to avoid Matlab recalculating the axes ticks
% you will need to set the following commands in you program:
%
% set(gca,'xtickmode','manual')
% set(gca,'ytickmode','manual')
% set(gca,'ztickmode','manual') % IF you have a third axis
%

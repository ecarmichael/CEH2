function Resize_figure(h, resize_factor)
%% Resize_figure: makes the figure bigger or smaller by the resize factor
%
%
%
%    Inputs: 
%     - h : figure handle.  default is gfc
%     - resize_factor: how much bigger (or smaller) default is 1.5
%
%
%
%    Outputs: 
%     - none
%
%
%
%
% EC 2020-01-15   initial version 
%
%
%% set default

if nargin == 0
    h = gcf;
    resize_factor = 1.5;
elseif nargin ==1
    resize_factor = 1.5;
end

if ~ishandle(h)
    error('Resize_figure: first input should be a figure handle.  try ''gcf''')
end
%% resize it
        pos = get(h, 'position');
    set(h,'position', [pos(1)-pos(1)/resize_factor, pos(2) pos(3)*resize_factor, pos(4)*resize_factor])
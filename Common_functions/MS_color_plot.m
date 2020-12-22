function h = MS_color_plot(tvec, dvec, type, cmap_in)
%% MS_color_plot: plot each data point in as a point in a color range. 
%
%
%
%    Inputs: 
%    - tvec  [1 x nSamples]    time vector
%
%    - dvec  [1 x nSamples]   data vector must be the same length as tvec
%
%    - type  [string]    type of plot markings '-', '.', 'x', ...
%
%    - cmap_in  [3 x n colors]   optional: an existing colormap. default is
%    parula.  
%
%    Outputs: 
%    - h     figure handle]
%
%
%
%
% EC 2020-12-20   initial version 
%
%
% % tried a piece of code from the great linspecer.m
% % Jonathan C. Lansey (2020). Beautiful and distinguishable line colors + colormap 
% % (https://www.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-colormap),
% % MATLAB Central File Exchange. Retrieved December 20, 2020.
%% initialize

if nargin <3
    type = '-'; 
    cmap_in = parula(length(tvec)); 
elseif nargin <4
    cmap_in = parula(length(tvec)); 
end


% %% set up color range (from linespecer.m)
%  x = linspace(1,length(tvec),size(cmap_in,1));
%     xi = 1:length(tvec);
%     cmap_out = zeros(length(tvec),3);
%     for ii=1:3
%         cmap_out(:,ii) = pchip(x,cmap_in(:,ii),xi);
%     end
% %     cmap_out = (cmap_out/255); 
%     
        
    %% plot each point
    hold on
    for ii = 1:length(tvec)
        plot(tvec(ii), dvec(ii), type, 'color', cmap_in(ii,:));
    end
    
    hold off
    
    h = gca; 
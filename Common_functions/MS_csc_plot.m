function h = MS_csc_plot(csc, offset)


if nargin < 2
    offset = 100; 

end




figure
clf

hold on

for ii = 1:length(csc.data)

    plot(csc.tvec, csc.data(ii,:)+ii*offset)
    y_t(ii) = mean(csc.data(ii,:)+ii*offset); 
    y_lab{ii} = csc.label{ii}; 
end


set(gca, 'YTick', y_t, 'YTickLabel', y_lab)


function h = sandbox_MS_Bar_scatter(d1, d2)
%% MS_Bar_scatter: makes a bar plot with the individual data points on top. 




%% test

d1 = randn(20,1); 
d2 = d1(1:end-2) + randn(1,1);



figure(909)
clf
bar(1:2, [mean(d1), mean(d2)])

hold on

    x1 =  1 + sort(MS_randn_range(length(d1), 1, -.2, .2)); 
    x2 = 2 + sort(MS_randn_range(length(d2), 1, -.2, .2)); 

    scatter(x1, d1, 25, 'r', 'filled')
    scatter(x2, d2, 25, 'b', 'filled')


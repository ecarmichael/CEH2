function yi = nearest_idx(x,y)
%% finds the y values with the minimal distance to x(ii).  based on nearest_idx2 from vandermeerlab github. 
yi = nan(size(x));

for ix = 1:length(x)
    
   d = abs(y-x(ix));
   [~,yi(ix)] = min(d);
    
end
function MS_Boris_NOL(fname)






%% init
disp(fname);


tbl = readtable(fname, 'VariableNamingRule', 'preserve');



%% check the file

behav = unique(tbl.Behavior); 
data = []; 

for ii = 1:length(behav)
    
    b_idx = contains(tbl.Behavior, behav{ii}); 
    
    data.types.(strrep(behav{ii}, ' ' , '_')) = tbl.("Behavior type")(b_idx); 
    types = unique( data.types.(strrep(behav{ii}, ' ' , '_'))); 
    
    
    fprintf('%s : ', behav{ii})
    for jj = 1:length(types)
        
        t_idx = contains(tbl.("Behavior type"), types{jj}) & b_idx;
        
        fprintf('%s = %2.0f evts |  ', types{jj}, sum(t_idx))
    end
    
    fprintf('\n')
    
end
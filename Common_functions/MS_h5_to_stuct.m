function out = MS_h5_to_stuct(fname)
%% get all the variable from an h5 file



%%  get the name of the datasets

info = h5info(fname);

out = []; 
for ii = 1:length(info.Datasets)
    fld = info.Datasets(ii).Name; 
    
   out.(fld) = h5read(fname, ['/' fld]); 
   
    
end
    


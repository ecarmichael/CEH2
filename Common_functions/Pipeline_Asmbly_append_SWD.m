function A_out = Pipeline_Asmbly_append_SWD(A_in, data_dir); 
%%  grab the LFP and get the SWR 

for iA = length(A_in):-1:1
   
load([data_dir filesep A_in{iA}.info.subject '_'  A_in{iA}.info.session '_LFP.mat'])


% for each assembly get the times and snap them to the LFP. 

A_in{iA}.REM_Pre_tvec
for ii = 1:length(this_LFP)
    
    
    REM_pre_csc.tvec = REM_pre_csc.tvec
    
end

end
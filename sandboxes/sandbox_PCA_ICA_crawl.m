%% sandbox_PCA_ICA_crawl

cd('C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter')

f_list = dir('pv*'); 

out = []; 
for ii = 4:length(f_list)
    
    [out.REM_z{ii}, out.Az{ii}] = sandbox_PCA_ICA(f_list(ii).name);
    
   clearvars -except f_list out ii 
   close all
end
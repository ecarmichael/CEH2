function A_in = Pipeline_Asmbly_append_preA(A_in)
%% check the pre_rem for assemblies, remove assemblies without positive weights, apply assemblies to the wake track data




%% get the pre REM templates
for ii = 1:length(A_in)
    
    if size(A_in{ii}.REM_Pre_data, 1) < size(A_in{ii}.REM_Pre_data, 2)
            fprintf('<strong>Pre REM data has less samples (%0.0f) than cells (%0.0f)</strong> \n',size(A_in{ii}.REM_Pre_data, 1), size(A_in{ii}.REM_Pre_data, 2))
            A_in{ii}.REM_temp = [];
            A_in{ii}.REM_Wake_proj = [];
        continue
    else
    
    REM_temp_all = assembly_patterns(A_in{ii}.REM_Pre_data); 
    
    Wake_proj = assembly_activity(REM_temp_all, A_in{ii}.wake_data'); 
    
    [REM_temp, Wake_proj]=  MS_Asmbly_select(REM_temp_all, Wake_proj, 2);
    
    A_in{ii}.REM_temp = REM_temp; 
    A_in{ii}.REM_Wake_proj = Wake_proj; 
    
    fprintf('[%.0f/%.0f = %.0f%%] Pre REM Assemblies had cells with positive weights (%.2fs binsize)\n',size(REM_temp,2),size(REM_temp_all,2),  (size(REM_temp,2)/size(REM_temp_all,2))*100, A_in{ii}.info.bin)
    end
end
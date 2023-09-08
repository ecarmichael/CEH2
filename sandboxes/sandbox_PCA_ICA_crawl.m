%% sandbox_PCA_ICA_crawl

cd('C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter')

f_list = dir('pv*');

out = []; 
for ii =1:length(f_list)

    [out.REM_z{ii}, out.Az{ii}] = sandbox_PCA_ICA(f_list(ii).name);
    
    clearvars -except f_list out ii
    close all
end

%%  same thing but stats only


% cd('C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter')
cd('/home/williamslab/Williams Lab Dropbox/Eric Carmichael/Comp_Can_inter')

f_list = dir('pv*');

out = [];
session = []; novel_idx = []; anx_idx = [];
for ii = 1:length(f_list)
    session{ii} = f_list(ii).name;
    [out.REM_z{ii}, out.Az{ii}, out.wake{ii}, out.rem{ii}] = sandbox_PCA_ICA_no_fig(f_list(ii).name);
    
    if ~isempty(strfind(f_list(ii).name, 'D1')) || ~isempty(strfind(f_list(ii).name, 'HATDS'))
        novel_idx(ii) = 1;
    end
    
    if isempty(strfind(f_list(ii).name, 'D5'))
        novel_idx(ii) = 0;
    end
    
    if ~isempty(strfind(f_list(ii).name, 'HAT'))
        anx_idx(ii) = 1;
    else
        anx_idx(ii) = 0;
    end
    
    
    %    clearvars -except f_list out ii
    %    close all
end

%% collect the output values
warning off
var_name = {'session', 'novel', 'anxiety', 'nAssemblies_z', 'nPlace_Assemblies', 'nWake_open', 'nWake_closed', 'nREM_open', 'nREM_closed', 'REM_OC_idx'};

for ii = length(session):-1:1
    nAssemblies_z = cell2mat(out.Az);
    
    if isstruct(out.REM_z{ii})
        nPlace_Assemblies(ii) = size(out.REM_z{ii}.all,1);
        nWake_open(ii) = sum(out.REM_z{ii}.isopen);
        nWake_close(ii) = sum(~out.REM_z{ii}.isopen);
        
        open_peaks = [];
        for jj = size(out.REM_z{ii}.open,1):-1:1
            [~, idx] = findpeaks(out.REM_z{ii}.open(jj,:), 'MinPeakHeight', 1, 'MinPeakDistance', 10);
            open_peaks(jj) = length(idx);
        end
        nREM_open(ii) = sum(open_peaks >0);
        
        
        
        % same thing for closed REM events
        
        closed_peaks = [];
        for jj = size(out.REM_z{ii}.close,1):-1:1
            [~, idx] = findpeaks(out.REM_z{ii}.close(jj,:), 'MinPeakHeight', 1, 'MinPeakDistance', 10);
            closed_peaks(jj) = length(idx);
        end
        nREM_close(ii) = sum(closed_peaks > 0);
        
        
        REM_OC_idx(ii) = (nREM_open(ii) - nREM_close(ii)) / size(out.REM_z{ii}.all,1);
    else
        nPlace_Assemblies = NaN;
        nWake_open(ii) = NaN;
        nWake_close(ii) = NaN;
        nREM_close(ii) = NaN;
        nREM_open(ii) = NaN;
        REM_OC_idx(ii) = NaN;
    end
end

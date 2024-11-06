function Master_SLEAP_NOL(data_dir)

cd(data_dir) 

a_list = dir('*.avi');
fnames = dir('label*analysis.h5');

if length(a_list) ~= length(fnames)
    
    error('# of .avi files and # of label h5 files are not equal')
end



v_num =[]; 
for ii = 1:length(fnames)

    p_idx = strfind(fnames(ii).name, '.')+1;
    n_idx(1) = strfind(fnames(ii).name(p_idx(1):p_idx(end)), '_')+p_idx(1);
    n_idx(2) = strfind(fnames(ii).name, '.analysis.h5');
    v_num(ii) = str2double(fnames(ii).name(n_idx(1):n_idx(2)-1));
    
    lab_list{ii} = fnames(ii).name;
end

[~, s_idx] = sort(v_num); 

lab_list_s = lab_list(s_idx);

%%

tsd = MS_SLEAP2TSD(lab_list_s, 30, [],[1, 1], 1)

%% plot to checl


plot(tsd.data(3,:), tsd.data(4,:), '.')

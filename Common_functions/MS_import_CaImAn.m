function ms = MS_import_CaImAn(data_dir, ts_name)

cnmfe_model = struct();


cnmfe_model.Groups(1).Datasets(6).Name = "C";
cnmfe_model.Groups(1).Datasets(6).Value = h5read(filename, "/estimates/C");
% 
% cnmfe_model.Groups(1).Datasets(14).Name = "S";
% cnmfe_model.Groups(1).Datasets(14).Value = h5read(fname, "/estimates/S");
% 
% 
% ms = []; 
% 
C = cnmfe_model.Groups(1).Datasets(6).Value; 
% ms.deconv = cnmfe_model.Groups(1).Datasets(14).Value; 

%% get the data using readNPY

ms = []; 

ms.RawTraces = readNPY([data_dir filesep 'dF_F.npy'])';
ms.deconv = readNPY([data_dir filesep 'C_calcium_traces.npy'])';
ms.denoise = readNPY([data_dir filesep 'S_spikes_deconvolved.npy'])';
SFPs = readNPY([data_dir filesep 'A_spatial_footprints.npy']);

ms.SFPs = [];
% convert the spf;
for ii = size(SFPs, 2):-1:1
    this_SFP = SFPs(:,ii);
    ms.SFPs(:,:,ii) = reshape(this_SFP, 600, 600);

    [~, idx] = max(ms.SFPs(:,:,ii), [], 'all'); 
    [ms.centroids(ii,2), ms.centroids(ii,1)] = ind2sub(size(ms.SFPs(:,:,ii)), idx);
end


ms = MS_append_timeStamps(ms, ts_name); 
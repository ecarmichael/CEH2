function ms = MS_append_deconv(ms)
%% MS_apped_deconv:  Applies the OASIS deconvulation method to ms data structure. Uses 'foopsi', 'ar2' (2nd order autoregressive)  

ms = msExtractBinary_detrendTraces(ms);
total_t = tic; 
if exist('deconvolveCa.m', 'file') == 2
    fprintf('\n<strong>%s</strong>: deconvolving traces...\n', mfilename)
    for iChan = size(ms.RawTraces,2):-1:1
            tic;
            [denoise,deconv] = deconvolveCa(ms.detrendRaw(:,iChan), 'foopsi', 'ar2', 'smin', -2.5, 'optimize_pars', true, 'optimize_b', true);
            toc;
            all_denoise(:,iChan) = denoise;    all_deconv(:,iChan) = deconv;
    end
    ms.denoise = all_denoise; 
    ms.deconv = all_deconv; 
    ms.decon_params = 'deconvolveCa(ms.detrendRaw(:,iChan), ''foopsi'', ''ar2'', ''smin'', -2.5, ''optimize_pars'', true, ''optimize_b'', true)'; 
    
%     total_toc = toc;
    toc(total_t);
    fprintf('<strong>%s</strong>: took %.2fs to process %.0f channels\n', mfilename, toc(total_t), size(ms.RawTraces,2));  
end
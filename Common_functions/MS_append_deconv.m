function ms = MS_append_deconv(ms, par_flag)
%% MS_apped_deconv:  Applies the OASIS deconvulation method to ms data structure. Uses 'foopsi', 'ar2' (2nd order autoregressive)


%% initialize

if nargin <2 
    par_flag = false;
end

%%

ms = msExtractBinary_detrendTraces(ms);
total_t = tic;
if ~exist('deconvolveCa.m', 'file') == 2
    error('OASIS deconvolveCa.m not found. Check that OASIS is in the path and oasis_setup has been run')
else
    fprintf('\n<strong>%s</strong>: deconvolving traces...\n', mfilename)
    
    all_denoise = NaN(size(ms.detrendRaw)); 
    all_deconv = all_denoise; 
    if par_flag
            detrendRaw = ms.detrendRaw; 

        parfor iChan = 1:size(ms.RawTraces,2)
            tic;
            [denoise,deconv] = deconvolveCa(detrendRaw(:,iChan), 'foopsi', 'ar2', 'smin', -2.5, 'optimize_pars', true, 'optimize_b', true);
            toc;
            all_denoise(:,iChan) = denoise;    all_deconv(:,iChan) = deconv;
        end
        
        
    else
        for iChan = size(ms.RawTraces,2):-1:1
            tic;
            [denoise,deconv] = deconvolveCa(ms.detrendRaw(:,iChan), 'foopsi', 'ar2', 'smin', -2.5, 'optimize_pars', true, 'optimize_b', true);
            toc;
            all_denoise(:,iChan) = denoise;    all_deconv(:,iChan) = deconv;
        end
    end
    ms.denoise = all_denoise;
    ms.deconv = all_deconv;
    ms.decon_params = 'deconvolveCa(ms.detrendRaw(:,iChan), ''foopsi'', ''ar2'', ''smin'', -2.5, ''optimize_pars'', true, ''optimize_b'', true)';
    
    %     total_toc = toc;
    toc(total_t);
    fprintf('<strong>%s</strong>: took %.2fs to process %.0f channels\n', mfilename, toc(total_t), size(ms.RawTraces,2));
end
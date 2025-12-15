function [ts, tsd] =  HF_piezo2ts(fname)
%% HF_piezo2iv: loads the Analog signal from the piezo, applies some filters and extracts the deflections as timestamps. 

        [data, tvec, info] = load_open_ephys_data(fname);
        csc = tsd(tvec, data);
        labels{ii} = info.header.channel;
        csc.cfg.hdr{ii} = info.header; 
        csc.cfg.hdr{ii}.SamplingFrequency = info.header.sampleRate; 
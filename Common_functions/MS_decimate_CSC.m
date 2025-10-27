function csc_out = MS_decimate_CSC(csc_in, factor)
%% :MS_decimate_CSC: decimate a TSD by factor.
%
%
%
%    Inputs:
%    - csc_in: [struct]   data in the TSD format
%    - factor: [double]   factor to decimate the data by. Should be a whole
%    number.
%
%
%    Outputs:
%    - csc_out: [struct]   decimated tsd struct.
%
%
%
%
% EC 2023-06-06   initial version
%
%
%
%% initialize


csc_out = csc_in;
csc_out.data = [];
csc_out.tvec = [];



fprintf('Decimating by factor %d old Fs:%d \n',factor, round(1/mode(diff(csc_in.tvec))))

for ii = size(csc_in.data,1):-1:1
    disp(csc_in.label{ii})
    csc_out.data(ii,:) = decimate(csc_in.data(ii,:),factor);
    csc_out.cfg.hdr{ii}.SamplingFrequency = csc_in.cfg.hdr{ii}.SamplingFrequency./factor;
    
end
csc_out.tvec = csc_in.tvec(1:factor:end);
fprintf('  new Fs: %d.\n', round(1/mode(diff(csc_out.tvec))))




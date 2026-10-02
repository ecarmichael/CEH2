function tsd_out = MS_select_tsd(tsd_in, idx);
%% MS_select_tsd: isolates channels of interest as a new tsd. 'idx' can be a number, index, or name;



% init

tsd_out = tsd_in;


% get the channels based on the idx type.

switch class(idx)

    case 'double'
        if isscalar(idx)
            channelIndex = idx;
        else
            channelIndex = idx(:);
        end

        ch_vec = ismember(1:size(tsd_in.data,1),idx);

    case 'cell'
        ch_vec = ismember(tsd_in.label,idx);
    case 'logical'
        if length(idx) == length(tsd_in.label)
            ch_vec
        else
            error('Logical index must match the number of channels.');
        end



    otherwise
    error('Unsupported index type: %s', class(idx));
end


%% remove the other channels

    tsd_out.data(~ch_vec,:) = [];
    tsd_out.label(~ch_vec) = [];

    tsd_out.cfg.hdr(~ch_vec) =[]; 

    CheckTSD(tsd_out)





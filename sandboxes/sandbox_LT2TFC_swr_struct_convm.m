%% sandbox_LT_tfc_inter_converter


all_sub = load('all_data_sub.mat');

all_sub = all_sub.all_data;



all_ca = load('all_data.mat');

all_ca = all_ca.all_data;

%% get common session

s_list = fieldnames(all_sub);
c_list = fieldnames(all_ca);

k_idx = ismember(c_list,s_list);

% s_list(~k_idx) = [];
c_list(~k_idx) = [];

% only keep the lt sessions
lt_idx= contains(s_list, 'tl');

s_list(~lt_idx) = [];
c_list(~lt_idx) = [];

%% loop and convert to individual files with both
% 2 and 7 don't have the right number of samples in the csc. 

for ii  = 1:length(s_list)
    if ii ==2  || ii == 7
        continue
    else
        data = [];

        data.(s_list{ii}).csc = all_ca.(s_list{ii}).csc;

        data.(s_list{ii}).csc.data(2,:) = all_sub.(s_list{ii}).csc.data;
        data.(s_list{ii}).csc.label{2} = all_sub.(s_list{ii}).csc.label{1};

        data.(s_list{ii}).swrs_ca1 = all_ca.(s_list{ii}).swrs;
        data.(s_list{ii}).swrs_sub = all_sub.(s_list{ii}).swrs;

        data.(s_list{ii}).evts = all_ca.(s_list{ii}).evts;

        if isfield(data.(s_list{ii}), 'S')
            data.(s_list{ii}).S = all_ca.(s_list{ii}).S;
        end

        if isfield(data.(s_list{ii}), 'ts_prime')
            data.(s_list{ii}).ts_prime = all_ca.(s_list{ii}).ts_prime;
        end

        % recompute the movement index
        move_ts = ts({all_ca.(s_list{ii}).evts.t{contains(all_ca.(s_list{ii}).evts.label, '8')}});

        mov_rate  = MS_spike2rate(move_ts, data.(s_list{ii}).csc.tvec);

        cfg_mov = [];
        cfg_mov.threshold = .001;
        cfg_mov.dcn = '<';
        cfg_mov.operation = '<';
        cfg_mov.minlen = .5;

        mov_iv = TSDtoIV(cfg_mov, mov_rate);

        % pad movement

        cfg_resize.d = [+.5 -.5];

        mov_iv = ResizeIV(cfg_resize, mov_iv);

        data.(s_list{ii}).mov_iv = mov_iv;
        data.(s_list{ii}).mov_vec = mov_rate;

        % save it back.
        save([strrep([s_list{ii}(1:3) s_list{ii}(5:end)], 'tl', 'LT') '.mat'], 'data')
    end

end
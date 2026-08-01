%% Pox_HF_LFP

% loops over sessions and get the the bandpower values for bands of
% interest. 


% load the Sub sessions. 

load('all_data_sub.mat')

s_list = fieldnames(all_data); 
pox_idx = nan(size(s_list)); 
pow_mat = NaN(length(s_list),5); 
psd = []; 

% index the pox mice(mouse)
for ii = 1:length(s_list)

    if sum(contains(s_list{ii}, '3265')) > 0
        pox_idx(ii) = ii+10;
    else
        pos_idx(ii) = 100;
    end
end

for ii = length(s_list):-1:1
   
    this_csc = all_data.(s_list{ii}).csc; 

    % get the movement from the encoder

    move_ts = ts({all_data.(s_list{ii}).evts.t{contains(all_data.(s_list{ii}).evts.label, '8')}});

    mov_rate  = MS_spike2rate(move_ts, this_csc.tvec);

    cfg_mov = [];
    cfg_mov.threshold = .001;
    cfg_mov.dcn = '<';
    cfg_mov.operation = '<';
    cfg_mov.minlen = .5;

    mov_iv = TSDtoIV(cfg_mov, mov_rate);

    % pad movement

    cfg_resize.d = [+.5 -.5];

    mov_iv = ResizeIV(cfg_resize, mov_iv);

    csc_mov = restrict(this_csc, mov_iv); 

    csc_no_mov = restrict(this_csc, InvertIV(mov_iv, this_csc.tvec(1), this_csc.tvec(end))); 


    %% grab the power from movement or non-movement epochs. 
    % overall
    pow_mat(ii,1) = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [1 40]);
    %theta
    pow_mat(ii,2) = bandpower(csc_mov.data, csc_mov.cfg.hdr{1}.SamplingFrequency, [6 10]); 
    % low gamma
    pow_mat(ii,3) = bandpower(csc_mov.data, csc_mov.cfg.hdr{1}.SamplingFrequency, [30 45]); 
    %mid gamma
    pow_mat(ii,4) = bandpower(csc_mov.data, csc_mov.cfg.hdr{1}.SamplingFrequency, [60 90]); 
    % ripples
    pow_mat(ii,5) = bandpower(csc_no_mov.data, csc_no_mov.cfg.hdr{1}.SamplingFrequency, [120 200]); 

    %% grab a psd for plotting. 

    hann_win = 2^12; % always make this in base 2 for speed

    [psd{ii}.pxx, psd{ii}.f] = pwelch(this_csc.data, hanning(hann_win), hann_win/2, hann_win*2 , this_csc.cfg.hdr{1}.SamplingFrequency);


end

%%  quick Psd comparisons

figure(101)
clf
hold on
ctrl_pxx = []; pox_pxx = []; 

for ii = length(psd):-1:1
    if pox_idx(ii)
    ctrl_pxx(ii,:) = 10*log10(psd{ii}.pxx)./max(10*log10(psd{ii}.pxx)); 
    else
    pox_pxx(ii,:) = 10*log10(psd{ii}.pxx)./max(10*log10(psd{ii}.pxx)); 

    end
end

pox_pxx(~pox_idx,:) = []; 
ctrl_pxx(pox_idx,:) = []; 

plot(psd{ii}.f, mean(ctrl_pxx))


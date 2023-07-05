%% sandbox SWR Starmaze


%% load the events and normalize


[evts_label, evts_ts, evts_hdr] = load_open_ephys_data('100_RhythmData-B.events'); 

evts_ts = evts_ts./30000; % normalize the original Fs. should always be 30k

% split out the different events



%% load data and decimate; 


cfg_load = [];
cfg_load.desired_sampling_frequency = 2000; 


csc = OE_old_csc2TSD(cfg_load);

%% normalize time to the start and stop of the Video

evts = OE_LoadEvents();


%% look at the signal for some time block
t1 = 0; t2 = 600;

t1_idx = nearest_idx3(t1, csc.tvec);
t2_idx = nearest_idx3(t2, csc.tvec);

c_ord = winter(length(csc.label)); 

figure(101)
clf
hold on


for ii  = length(csc.label):-1:1
    
    plot(csc.tvec(t1_idx:t2_idx), csc.data(ii, t1_idx:t2_idx) +ii*1000, 'color', c_ord(ii,:)); 
    tick_y(ii) = median(csc.data(ii, t1_idx:t2_idx) +ii*1000); 
    ch_val = regexp(csc.label{ii},'\d*','Match'); 
    ch_num(ii) = str2num(ch_val{1});
end

xlim([csc.tvec(t1_idx), csc.tvec(t2_idx)]); 

set(gca, 'YTick', tick_y, 'YTickLabel', csc.label); 

y_lims = ylim; 
%% append events to the trace plot

lin_ord = linspecer(length(evts.t)); 

for iE = 2:length(evts.t)
   
%     plot([evts.t{iE}, evts.t{iE}]', repmat(y_lims, length(evts.t{iE}), 1)','color', lin_ord(iE,:))
    
    for ii  = length(evts.t{iE}):-1:1
        h = vline(evts.t{iE}(ii), 'k', evts.label{iE}); 
        h.Color = lin_ord(iE,:); 
    end
    
    
end

%% trim to TTL onset
TTL_idx = find(contains(evts.label, '6')); 

start_ts = evts.t{TTL_idx}(1); 
end_ts = evts.t{TTL_idx}(end); 


csc_r = restrict(csc, start_ts, end_ts); 

csc_r.tvec = csc_r.tvec - csc_r.tvec(1); 

evts_r = restrict(evts, start_ts, end_ts); 

for iE = 1:length(evts_r.t)
    evts_r.t{iE} = evts_r.t{iE} - start_ts; 
end

%% good channels
good_SWR = [3, 21, 32, 33, 35, 59, 63]; 
keep_idx = ismember(ch_num, good_SWR);

csc.data(~keep_idx, :) = [];
csc.label(~keep_idx) = []; 
csc.cfg.hdr(~keep_idx) = []; 

%% filter down to 3 - 400 hz
cfg_filt.check = 0; % plot checks.
cfg_filt.filt.type = 'butter'; %Cheby1 is sharper than butter
cfg_filt.filt.f  = [4 400]; % broad, could use 150-200?
cfg_filt.filt.order = 4; %type filter order (fine for this f range)
cfg_filt.filt.display_filter = 0; % use this to see the fvtool

csc_filt = FilterLFP(cfg_filt, csc); 


t1 = 0; t2 = 600;

t1_idx = nearest_idx3(t1, csc.tvec);
t2_idx = nearest_idx3(t2, csc.tvec);

c_ord = winter(length(csc.label)); 

figure(101)
clf
hold on

for ii  = length(csc.label):-1:1
    
    plot(csc_filt.tvec(t1_idx:t2_idx), csc_filt.data(ii, t1_idx:t2_idx) +ii*1000, 'color', c_ord(ii,:)); 
    tick_y(ii) = median(csc_filt.data(ii, t1_idx:t2_idx) +ii*1000); 
    ch_val = regexp(csc_filt.label{ii},'\d*','Match'); 
    ch_num(ii) = str2num(ch_val{1});
end

xlim([csc_filt.tvec(t1_idx), csc_filt.tvec(t2_idx)]); 
set(gca, 'YTick', tick_y, 'YTickLabel', csc.label); 

%% run the SWR on a few channels and then look for intersects. 

close all
SWR.('CH33') = MS_SWR_detector(csc, 'CH33',0)

pause(2)
close all

SWR.('CH21') = MS_SWR_detector(csc, 'CH21',0)

pause(2)
close all

%% look for intersecting SWRs across a pair of channels
cfg_int = [];
% cfg_int.
SWR_evts = MS_IntersectIV(cfg_int, SWR.CH21, SWR.CH33);


close all
cfg_plot= [];
    cfg_plot.display = 'iv'; %'iv';
    cfg_plot.title = 'num';
    cfg_plot.target = 'CH33'; 
PlotTSDfromIV(cfg_plot, SWR_evts, csc_filt)


%% remove bad events
% Baseline
% rm_idx = [4 8 11 25 30 32  152 157 158 161 166 167 169 172 173 176 177 178 179 180 181 183 184 186 187 188 189  190 191 194 198 206 207 208 209 211 213 214 217 218 219 223 226 237 238 240 242 243 244 266 267 269 274 275 277 278 279 280 281 282 284 291 296 297 298  299 300 301 303 307 308 309 311 312 313 317 320 322 324 325 326 328 334 336 338 339 340 347 348 349 353 354 355 357 358 366 370 375 400 401 402 403 406];

% Omega 1
rm_idx = [37 47 51 58 88 100 103 104 113 132 135 163 169 170 173 202 213 223 237 251 252 266 282 322 333 385 408 417 424 474 476 482 504 516 519 521 523 532 535 551 552 563 566 569 577 579 603 632 633]; 

rm_idx = ismember(SWR_evts.usr.num, rm_idx); 

SWR_out = SelectIV([], SWR_evts, ~rm_idx); 


%% final check
close all
cfg_plot= [];
    cfg_plot.display = 'iv'; %'iv';
    cfg_plot.title = 'num';
    cfg_plot.target = 'CH35'; 
PlotTSDfromIV(cfg_plot, SWR_out, csc_filt)


%     swr_cent = IVcenters(SWR_out);
    win = csc.cfg.hdr{1}.SamplingFrequency / 5;
    
    swr_idx = nearest_idx3(SWR_out.tstart, csc.tvec);
    
    all_SWR = nan(length(SWR_out.tstart),win*2+1) ;
    for ii =length(swr_cent):-1:1
        if (swr_idx(ii) + win) < length(csc.tvec)
            all_SWR(ii,:) = csc_filt.data(1,swr_idx(ii) - win: swr_idx(ii)+win);
        end
    end
    
    
    figure(909)
    
    subplot(2,1,1)
    
    plot(-.2:1/(csc.cfg.hdr{1}.SamplingFrequency):.2, nanmean(all_SWR))
    
    
    figure
    subplot(2,1,2)
    csc_chan = contains(csc_filt.label, cfg_plot.target); 
    csc_spec = csc_filt; 
    csc_spec.data(~csc_chan,:) = []; 
    csc_spec.label(~csc_chan) = []; 
    csc_spec.cfg.hdr(~csc_chan) = []; 
    
    Triggered_Spec_FT(csc_spec, IVcenters(SWR_out), 'SWR', 100:.1:250, [],  [-.1 .1])




%% split out events that are singles, doubles, or triples. 

% from Yamamoto & Tonegawa 2017 (https://www.sciencedirect.com/science/article/pii/S0896627317308577#sec4)
% "After detecting candidate ripple events (see Candidate Ripple
% and Replay Event Detection), time lags between peaks of ripple event were computed. 
% Singlet ripples were defined as ripple events that are temporally separated by 200 ms 
% or more to adjacent ripple events. For the doublets and triplets (ripple bursts), 
% the temporal lags were set to less than 200 ms but greater than 70 ms. 
% These parameters were determined by computing the single ripple duration 
% (∼80 ms) and inter-ripple intervals (∼125 ms) when ripple bursts occurred. 
% In the event that adjacent ripple closer than 70 ms, 
% these were categorized as single ripple. Number of adjacent ripple peaks was 
% determined and assigned to doublets, triplets and others (more peaks than
% three). "

% cfg_merge = [];
% cfg_merge.gap = .07; 
% SWR_merge = MergeIV(cfg_merge, SWR_out); 


SWR_IPI = diff(IVcenters(SWR_out)); 

merge_idx = SWR_IPI < 0.07; 





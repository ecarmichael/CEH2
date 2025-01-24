function  out = MS_DLC_score_freezing(fname, proto, LED_on)
%% MS_score_freezing: score freezing based on movement in DLC tracking data. if multiple body parts, average across to get a better measure. 









%% initialize 
if nargin < 2
    proto = []; 
end


% find the label file
f_list = dir('*filtered.csv');

keep_idx = zeros(1, length(f_list)); 

for f = 1:length(f_list)
    if contains(f_list(f).name, fname(1:end-4))
        keep_idx(f) = true; 
    else
        keep_idx(f) = false;
    end
end

keep_idx = logical(keep_idx); 

dlc_name = f_list(keep_idx).name; 
dlc_dir = f_list(keep_idx)

if contains(fname, 'TFC2')
    conv_fact = [47.5, 44];  %
else
    conv_fact = [47.5, 44]; % convert from pixels (890 x 800) to cm (box dimensions [48 44])
end

pos = MS_DLC2TSD_single(dlc_name, fname,conv_fact); 

%trim to the LED off signal. 
s_idx = int16(LED_on+(2/mode(diff(pos.tvec)))); 
e_idx = int16(LED_on+(2/mode(diff(pos.tvec))) + proto.ITI(end,2)/mode(diff(pos.tvec)));
pos_r = pos; 
pos_r.tvec  = pos.tvec(s_idx:e_idx); 
pos_r.tvec = pos_r.tvec - pos_r.tvec(1); 
pos_r.data  = pos.data(:,s_idx:e_idx); 

%% extract the movement from the head speed vector. 

% get the moving mean of the speed
t_win = 2; 
t_win_f = round(t_win /mode(diff(pos_r.tvec))); % convert to frames/samples

mov_m = movmedian(pos_r.data(end-1,:), t_win_f); 

fvec = zeros(1, length(mov_m)); 

fvec(mov_m < prctile(mov_m, 50)) = 1; 

fvec = logical(fvec); 

pos_f = pos_r; 
pos_f.data = [];
pos_f.label = []; 
pos_f.data(end+1,:) = fvec; 
pos_f.label{end+1} = 'freeze_vec';

% remove blocks with less than 2s freezing. 
cfg_f = []; 
cfg_f.method = 'raw';
cfg_f.threshold = .5;
cfg_f.operation = '>'; 
cfg_f.minlen = 2; 

f_iv = TSDtoIV(cfg_f, pos_f); 


figure(101)
clf
subplot(1,4,1:3)

hold on

plot(pos_r.tvec, pos_r.data(end-1,:), 'k')
plot(pos_r.tvec, mov_m, 'b')

yline(prctile(mov_m, 50), 'r', '50th prctile move median')

plot(pos_r.tvec,fvec)

proto_fnames = fieldnames(proto);
c_ord = MS_linspecer(length(proto_fnames)); 


%%
for F = 1:length(fieldnames(proto))
    
    evts = size(proto.(proto_fnames{F}),1);
    
    for ii = 1:evts; 
        
        rectangle('position', [proto.(proto_fnames{F})(ii,1), -5, proto.(proto_fnames{F})(ii,2) - proto.(proto_fnames{F})(ii,1), 5], 'FaceColor', c_ord(F,:))
        text((proto.(proto_fnames{F})(ii,1) +( proto.(proto_fnames{F})(ii,2) - proto.(proto_fnames{F})(ii,1))/2), -2.5, proto_fnames{F}, 'HorizontalAlignment', 'center')
    end
end

%% get the minute by minue
f_val = zeros(1, length(0:2:ceil(pos_r.tvec(end))); 

for ii = 0:2:ceil(pos_r.tvec(end))
    
    if ii == 0
        s_idx = 1; 
    else
        
    s_idx = int16((ii*360)/mode(diff(pos_r.tvec))); 
    end
    
    if ii == ceil(pos_r.tvec(end) / 60)-1
        e_idx = length(pos_r.tvec);
    else
            e_idx = int16(((ii+1)*360)/mode(diff(pos_r.tvec)));
    end
    
    f_val(ii+1) = sum(fvec(s_idx:e_idx))/length(fvec(s_idx:e_idx)); 
    fprintf('block length %f  F: %.2f %%\n', (e_idx - s_idx)*mode(diff(pos_r.tvec)),f_val(ii+1))
    
end

figure(101)
subplot(1,4,4)
plot(0:2:ceil(pos_r.tvec(end))+1, f_val)





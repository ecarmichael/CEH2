function  out = MS_DLC_score_freezing_dir(fname, thresh, proto, LED_on, fig_dir)
%% MS_score_freezing: score freezing based on movement in DLC tracking data. if multiple body parts, average across to get a better measure.









%% initialize
if nargin < 2
    thresh = 50;
    proto = [];
    LED_on = 1;
    fig_dir = [];
elseif nargin < 3
    proto = [];
    LED_on = 1;
    fig_dir = [];
elseif nargin < 4
    LED_on = 1;
    fig_dir = [];
elseif nargin < 5
    fig_dir = [];
    
end


% find the label file
f_list = dir('*filtered.csv');

if contains(fname, 'TFC2')
    conv_fact = [1125/46, 668/29];  % 46 x 29cm. (1125 x 668pixels)
else
    conv_fact = [47.5, 44]; % convert from pixels (890 x 800) to cm (box dimensions [48 44])
end

pos = MS_DLC2TSD(fname,[], conv_fact);

%trim to the LED off signal.
if contains(fname, 'TFC3') || contains(fname, 'TFC_C')
    s_idx = nearest_idx(0, pos.tvec);
    e_idx = nearest_idx(proto.baseline(end), pos.tvec);
else
    s_idx = int16(LED_on);
    e_idx = s_idx +nearest_idx(proto.ITI(end,2), pos.tvec);
end

pos_r = restrict(pos, pos.tvec(s_idx), pos.tvec(end));
pos_r.tvec  = pos.tvec(s_idx:end);
pos_r.tvec = pos_r.tvec - pos_r.tvec(1);
pos_r.data  = pos.data(:,s_idx:end);

% vx = dxdt(pos_r.tvec,pos_r.data(end-3,:));
% vy = dxdt(pos_r.tvec,pos_r.data(end-2,:));
% 
% pos_r.data(end+1,:) = sqrt(vx.^2+vy.^2)';
% pos_r.label{end+1} = 'speed body'; 
% 
% spd_idx = find(contains(pos_r.label, 'Speed'));

spd_idx = size(pos_r.data,1) -1;
%% extract the movement from the head speed vector.

% get the moving mean of the speed
t_win = 2.5;
t_win_f = round(t_win /mode(diff(pos_r.tvec))); % convert to frames/samples

mov_m = movmedian(pos_r.data(spd_idx,:), t_win_f);

fvec = zeros(1, length(mov_m));

if isempty(thresh)
    mov_thresh = .5;
else
    mov_thresh = thresh; 
%     mov_thresh = prctile(mov_m, thresh);
end

fvec(mov_m < mov_thresh) = 1;

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

if ~isempty(fig_dir)
figure(101)
clf
ax(1) = subplot(3,2,[1:4]);

hold on

plot(pos_r.tvec, pos_r.data(spd_idx,:), 'k')
plot(pos_r.tvec, mov_m, 'b')

if isempty(thresh)
    yline(mov_thresh, 'r', [num2str(thresh) 'th prctile move median'])
    
else
    yline(prctile(mov_m, thresh), 'r', [num2str(mov_thresh) 'th prctile move median'])
end
plot(pos_r.tvec,fvec)

proto_fnames = fieldnames(proto);
c_ord = MS_linspecer(length(proto_fnames));

end
%%
for F = 1:length(fieldnames(proto))
    
    evts = size(proto.(proto_fnames{F}),1);
    this_mean = [];
    
    for ii = 1:evts
        
        out.TFC.(proto_fnames{F})(ii) = sum(fvec(nearest_idx(proto.(proto_fnames{F})(ii,1), pos_r.tvec):nearest_idx(proto.(proto_fnames{F})(ii,2), pos_r.tvec)))...
            /length(nearest_idx(proto.(proto_fnames{F})(ii,1), pos_r.tvec):nearest_idx(proto.(proto_fnames{F})(ii,2), pos_r.tvec));
        
        if ~isempty(fig_dir)

        rectangle('position', [proto.(proto_fnames{F})(ii,1), -5, proto.(proto_fnames{F})(ii,2) - proto.(proto_fnames{F})(ii,1), 5], 'FaceColor', c_ord(F,:))
        text((proto.(proto_fnames{F})(ii,1) +( proto.(proto_fnames{F})(ii,2) - proto.(proto_fnames{F})(ii,1))/2), -2.5, [proto_fnames{F} ' ' num2str(out.TFC.(proto_fnames{F})(ii)*100,3) '%'], 'HorizontalAlignment', 'center')
        end
        out.TFC.(proto_fnames{F})(ii) = sum(fvec(nearest_idx(proto.(proto_fnames{F})(ii,1), pos_r.tvec):nearest_idx(proto.(proto_fnames{F})(ii,2), pos_r.tvec)))...
            /length(nearest_idx(proto.(proto_fnames{F})(ii,1), pos_r.tvec):nearest_idx(proto.(proto_fnames{F})(ii,2), pos_r.tvec));
    end
end

%% get the minute by minue

t_bin =10;

if contains(fname, 'TFC3') || contains(fname, 'TFC_C')
    cut = 300; 
else
    cut = proto.ITI(end,2); 

end
    f_val = zeros(1, length(0:t_bin:cut-t_bin));

c= 0;
for ii = 0:t_bin:cut-t_bin
    c = c+1;
    
    s_idx = nearest_idx(ii, pos_r.tvec);
    
    %     if ii == ceil(pos_r.tvec(end))-2
    %         e_idx = length(pos_r.tvec);
    %     else
    e_idx = nearest_idx(ii+t_bin,pos_r.tvec);
    %     end
    
    f_val(c) = sum(fvec(s_idx:e_idx))/length(fvec(s_idx:e_idx));
    %fprintf('block length %.0f  F: %.2f%%\n', (e_idx - s_idx)*mode(diff(pos_r.tvec)),f_val(c)*100)
    
end

if ~isempty(fig_dir)

figure(101)
ax(2) = subplot(3,2,[5:6]);
bar((0:t_bin:cut-t_bin)+(t_bin/2), f_val,'BarWidth', 1)
% set(gca, 'XTicklabel', (0:t_bin:ceil(pos_r.tvec(end))-2))


linkaxes(ax, 'x')
xlim([pos_r.tvec(1) pos_r.tvec(end)]);

saveas(gcf, fig_dir)
end

%% collect the outputs
out.pos_r = pos_r; 
out.fvec = fvec;
out.f_bin = f_val;
out.t_bin = (0:t_bin:cut-t_bin);


function [emp_idx, labels, boxes] = MS_DLC_EMP(data_dir, save_dir, boxes, plt_mov)

if nargin < 2
    save_dir = [];
    boxes = [];
    plt_mov = 0;
elseif nargin <3
    boxes = [];
    plt_mov = 0;
elseif nargin <4
    plt_mov = 0;
    
end

cd(data_dir)

a_list = dir('*.avi');
fnames = dir('*filtered.csv');

if length(a_list) ~= length(fnames)
    
    error('# of .avi files and # of label h5 files are not equal')
end



%% Load and scale the position data from DLC

conv_f = [1 1];

tsd = MS_DLC2TSD(data_dir, [],conv_f, 0);

if (tsd.tvec(end) - tsd.tvec(1)) > 600
    
    s_idx = nearest_idx(5, tsd.tvec);
    e_idx = nearest_idx(605, tsd.tvec);
    
    tsd.tvec = tsd.tvec(s_idx:e_idx);
    tsd.data = tsd.data(:,s_idx:e_idx);
    
    
end

%% plot to checl




%% set the boarders and collect the
c_ord = MS_linspecer(5);

figure(191)
clf
hold on

imagesc(tsd.mean_frame)
plot(tsd.data(3,:), tsd.data(4,:), '.');

if isempty(boxes)
    
    closed = MS_drawrectangle_wait('Closed',c_ord(1,:));
    open_trk = MS_drawrectangle_wait('Open track',[1 1 0]);
    open_wide_R = MS_drawrectangle_wait('Open wide right',[0 1 1]);
    open_wide_L = MS_drawrectangle_wait('Open wide left',[1 0 1]);
    
    [~, ~, close_p] = MS_rec2corner(closed.Position);
    [~, ~, open_p] = MS_rec2corner(open_trk.Position);
    [~, ~, open_wR_p] = MS_rec2corner(open_wide_R.Position);
    [~, ~, open_wL_p] = MS_rec2corner(open_wide_L.Position);
    
else
    close_p = boxes.close;
    open_p = boxes.open;
    open_wR_p = boxes.open_wR;
    open_wL_p = boxes.open_wL;
end


xcross = intersect(open_p,close_p);

hold on
plot(xcross, 'FaceColor', [0 1 0], 'FaceAlpha', .5, 'EdgeColor', [0 1 0], 'LineWidth', 5);

%% vectorize occupancy

head_idx = [3 4];
body_idx = [5 6]; %Which body node to use.

% find points where the head is not close to the body

e_dist = sqrt(((tsd.data(head_idx(1),:) - tsd.data(body_idx(1),:)).^2)+((tsd.data(head_idx(2),:) - tsd.data(body_idx(2),:)).^2));


%
% figure(878)
% clf
% imagesc(tsd.mean_frame)
% colormap('bone')
%
% hold on
%
% plot(close_p, 'FaceColor',[ 0 0 1], 'EdgeColor',[0 0 1], 'FaceAlpha', .2)
% plot(open_p, 'FaceColor',[ 1 1 0], 'EdgeColor',[1 1 0], 'FaceAlpha', .2)
% plot(xcross, 'FaceColor', [0 1 0], 'FaceAlpha', .5, 'EdgeColor', [0 1 0], 'LineWidth', 3);
% plot(open_wR_p, 'FaceColor', [0 1 1], 'FaceAlpha', .5, 'EdgeColor', [0 1 1], 'LineWidth', 3);
% plot(open_wL_p, 'FaceColor', [1 0 1], 'FaceAlpha', .5, 'EdgeColor', [1 0 1], 'LineWidth', 3);


[c_in, c_on] = inpolygon(tsd.data(body_idx(1),:), tsd.data(body_idx(2),:), close_p.Vertices(:,1), close_p.Vertices(:,2));
c_idx = c_in | c_on;

% [oWR_in, oWR_on] = inpolygon(tsd.data(body_idx(1),:), tsd.data(body_idx(2),:), open_wR_p.Vertices(:,1), open_wR_p.Vertices(:,2));
%
% [oWL_in, oWL_on] = inpolygon(tsd.data(body_idx(1),:), tsd.data(body_idx(2),:), open_wL_p.Vertices(:,1), open_wL_p.Vertices(:,2));

[o_in, o_on] = inpolygon(tsd.data(body_idx(1),:), tsd.data(body_idx(2),:), open_p.Vertices(:,1), open_p.Vertices(:,2));
% o_idx = o_in | o_on | oWR_in | oWR_on | oWL_in | oWL_on;
o_idx = o_in | o_on;


[x_in, x_on] = inpolygon(tsd.data(body_idx(1),:), tsd.data(body_idx(2),:), xcross.Vertices(:,1), xcross.Vertices(:,2));
x_idx = x_in | x_on;

c_idx = c_idx & (~x_idx & ~o_idx);
o_idx = o_idx & (~x_idx & ~c_idx);
n_idx = ~x_idx & ~o_idx & ~c_idx;

% get putative head dips.
[h_o_in, h_o_on] = inpolygon(tsd.data(head_idx(1),:), tsd.data(head_idx(2),:), open_p.Vertices(:,1), open_p.Vertices(:,2));

[h_wr_in, h_wr_on] = inpolygon(tsd.data(head_idx(1),:), tsd.data(head_idx(2),:), open_wR_p.Vertices(:,1), open_wR_p.Vertices(:,2));
[h_wl_in, h_wl_on] = inpolygon(tsd.data(head_idx(1),:), tsd.data(head_idx(2),:), open_wL_p.Vertices(:,1), open_wL_p.Vertices(:,2));

d_idx = (h_wr_in | h_wr_on | h_wl_in | h_wl_on) & ~(h_o_in | h_o_on)  & ~(c_idx | x_idx);


c_idx = c_idx & ~d_idx;


% p1 = plot(tsd.data(body_idx(1),c_idx),tsd.data(body_idx(2),c_idx),'b+'); % points inside
% p2 = plot(tsd.data(body_idx(1),o_idx),tsd.data(body_idx(2),o_idx),'yo'); % points outside
% p3 = plot(tsd.data(body_idx(1),x_idx),tsd.data(body_idx(2),x_idx),'gs'); % points outside
% p4 = plot(tsd.data(body_idx(1),n_idx),tsd.data(body_idx(2),n_idx),'kd'); % points outside
% p5 = plot(tsd.data(head_idx(1),d_idx),tsd.data(head_idx(2),d_idx),'rx'); % points outside





% lg.String(1:5)  = [];
%% collect the data

emp_idx = NaN(size(tsd.data(3,:)));

emp_idx(n_idx) = 0;
emp_idx(c_idx) = 1;
emp_idx(o_idx) = 2;
emp_idx(x_idx) = 3;

% [c_entries] = find(diff(emp_idx == 1) > 0);
% [o_entries] = find(diff(emp_idx == 2) > 0);
% [x_entries] = find(diff(emp_idx == 3) > 0);
%
% [c_exits] = find(diff(emp_idx == 1) < 0);
% [o_exits] = find(diff(emp_idx == 2) < 0);
% [x_exits] = find(diff(emp_idx == 3) < 0);

%% get transitons
% %%%%% Closed to Transition %%%%%%%%%
c_entries = find(diff(c_idx == 1) >0);
c_exits    = find(diff(c_idx == 1) <0);

if c_entries(1) > c_exits(1)
    c_entries = [1 c_entries]; 
end
% make 'trials' for each entry/exit type.
if length(c_entries) > length(c_exits)
    c_exits(end+1) = length(c_idx);
end

% use IV functions to link transitions or something.
CtT_IV = iv(c_entries, c_exits);

cfg_m = [];
cfg_m.gap = .5* round(1/mode(diff(tsd.tvec)));
CtT_IV_m = MergeIV(cfg_m, CtT_IV);


% %%%%%transition 2 open %%%%%%%%%
o_entries = find(diff(o_idx == 1) >0);
o_exits    = find(diff(o_idx == 1) <0);

if o_entries(1) > o_exits(1)
    o_entries = [1 o_entries]; 
end
% make 'trials' for each entry/exit type.
if length(o_entries) > length(o_exits)
    o_exits(end+1) = length(o_idx);
end

% use IV functions to link transitions or something.
OtT_IV = iv(o_entries, o_exits);
cfg_m = [];
cfg_m.gap = .5* round(1/mode(diff(tsd.tvec)));
OtT_IV_m = MergeIV(cfg_m, OtT_IV);


% %%%%% get head dips %%%%%%%%%%%%
d_entries = find(diff(d_idx == 1) >0);
d_exits    = find(diff(d_idx == 1) <0);
if d_entries(1) > d_exits(1)
    d_entries = [1 d_entries]; 
end
% make 'trials' for each entry/exit type.
if length(d_entries) > length(d_exits)
    d_exits(end+1) = length(d_idx);
end

% use IV functions to link transitions or something.
dip_IV = iv(d_entries, d_exits);

cfg_m = [];
cfg_m.gap = .5* round(1/mode(diff(tsd.tvec)));
dip_IV_m = MergeIV(cfg_m, dip_IV);

d_entries = dip_IV_m.tstart;
d_exits = dip_IV_m.tend;

d_idx = zeros(size(n_idx));
for ii = 1:length(d_entries)
    d_idx(d_entries(ii):d_exits(ii)) = 1;
end
d_idx = logical(d_idx);

emp_idx(d_idx) = 4;


open_prct = (sum(o_idx)/length(o_idx))*100;
closed_prct = (sum(c_idx)/length(c_idx))*100;
xcross_prct = (sum(x_idx)/length(x_idx))*100;
h_dips_prct = (sum(d_idx)/length(d_idx))*100;
none_prct = (sum(n_idx)/length(n_idx))*100;

fprintf('EMP: %0.1f%% open | %0.1f%% closed | %0.1f%% xcross | %0.1f%% head dips| %0.1f%% none\n',open_prct,  closed_prct, xcross_prct,h_dips_prct, none_prct)

figure(8818)
clf
subplot(2,2,2)
p = pie([(sum(c_idx)/length(c_idx)),(sum(o_idx)/length(o_idx)), (sum(x_idx)/length(x_idx)),(sum(d_idx)/length(d_idx))],[0 1 1 1],'%.3f%%');

c = 0;
for ip = 1:2:length(p)
    c = c+1;
    p(ip).FaceColor = c_ord(c,:);
    
end

legend({'Open' , 'Closed', 'xCross', 'Head dip'}, 'location', 'eastoutside')

title(' Percentage of time in each region')

subplot(2,2,3:4)
hold on
plot(tsd.tvec,emp_idx-0.5, 'k', 'linewidth', .5);
plot(tsd.tvec(c_idx), emp_idx(c_idx)-.5,'.', 'color', c_ord(1,:), 'markersize', 22)
plot(tsd.tvec(o_idx), emp_idx(o_idx)-.5,'.', 'color', c_ord(2,:), 'markersize', 22)
plot(tsd.tvec(x_idx), emp_idx(x_idx)-.5,'.', 'color', c_ord(3,:), 'markersize', 22)
plot(tsd.tvec(n_idx), emp_idx(n_idx)-.5,'.', 'color', c_ord(4,:), 'markersize', 22)
plot(tsd.tvec(d_idx), emp_idx(d_idx)-.5,'.', 'color', c_ord(5,:), 'markersize', 22)
xlim([tsd.tvec(1), tsd.tvec(end)])
set(gca,'ytick', -.5:3.5, 'yticklabel', {'none', 'closed', 'open', 'cross', 'head dip'})
hline(0:3);
ylim([-1 4])
% transitions
min_t = 2; % time in seconds for the animal to occuppy a region for it to count.
min_t = min_t * round(1/mode(diff(tsd.tvec)));

subplot(2,2,1)
hold on

imagesc(tsd.mean_frame)
colormap('bone')
axis off
hold on

plot(close_p, 'FaceColor',c_ord(1,:), 'EdgeColor',c_ord(1,:), 'FaceAlpha', .1)
plot(open_p, 'FaceColor',c_ord(2,:), 'EdgeColor',c_ord(2,:), 'FaceAlpha', .1)
plot(xcross, 'FaceColor', c_ord(3,:), 'FaceAlpha', .1, 'EdgeColor', c_ord(3,:), 'LineWidth', 3);
plot(open_wR_p, 'FaceColor', c_ord(5,:), 'FaceAlpha', .1, 'EdgeColor', c_ord(5,:), 'LineWidth', 3);
plot(open_wL_p, 'FaceColor', c_ord(5,:), 'FaceAlpha', .1, 'EdgeColor', c_ord(5,:), 'LineWidth', 3);

p1 = plot(tsd.data(body_idx(1),c_idx),tsd.data(body_idx(2),c_idx),'.','color', c_ord(1,:), 'markersize', 12); % points inside
p2 = plot(tsd.data(body_idx(1),o_idx),tsd.data(body_idx(2),o_idx),'o','color',  c_ord(2,:), 'markersize', 12); % points outside
p3 = plot(tsd.data(body_idx(1),x_idx),tsd.data(body_idx(2),x_idx),'x','color',  c_ord(3,:), 'markersize', 12); % points outside
p4 = plot(tsd.data(body_idx(1),n_idx),tsd.data(body_idx(2),n_idx),'.','color',  c_ord(4,:), 'markersize', 12); % points outside
p5 = plot(tsd.data(head_idx(1),d_idx),tsd.data(head_idx(2),d_idx),'s','color',  c_ord(5,:), 'markersize', 12); % points outside

lg = legend([p1 p2 p3 p4 p5], {'closed idx', 'open idx', 'x idx', 'none idx', 'dip idx'});


labels  = {'none', 'closed', 'open', 'trans', 'head dip'};
axis equal

%% debug movie

if plt_mov
    figure(919)
    clf
    imagesc(tsd.mean_frame)
    colormap('bone')
    axis off
    hold on
    
    plot(close_p, 'FaceColor',c_ord(1,:), 'EdgeColor',c_ord(1,:), 'FaceAlpha', .1)
    plot(open_p, 'FaceColor',[ 1 1 0], 'EdgeColor',[1 1 0], 'FaceAlpha', .1)
    plot(xcross, 'FaceColor', [0 1 0], 'FaceAlpha', .1, 'EdgeColor', [0 1 0], 'LineWidth', 3);
    plot(open_wR_p, 'FaceColor', [0 1 1], 'FaceAlpha', .1, 'EdgeColor', [0 1 1], 'LineWidth', 3);
    plot(open_wL_p, 'FaceColor', [1 0 1], 'FaceAlpha', .1, 'EdgeColor', [1 0 1], 'LineWidth', 3);
    
    
    m_t = [min(tsd.data(1,:)) min(tsd.data(2,:))];
    % mx_t =
    
    text(m_t(1), m_t(2), min(tsd.data(2,:))+((max(tsd.data(2,:)) - min(tsd.data(2,:)))/5),['10x speed'], 'color', 'w')
    
    
    for iF =length(emp_idx):-1:1
        
        
        p1 = plot(tsd.data(body_idx(1),iF),tsd.data(body_idx(2),iF),'s','color', c_ord(emp_idx(iF),:), 'markersize', 12); %
        p2 = plot(tsd.data(head_idx(1),iF),tsd.data(head_idx(2),iF),'o','color', c_ord(emp_idx(iF),:), 'markersize', 12); %
        
        t1  = text(min(tsd.data(1,:)), min(tsd.data(2,:))+((max(tsd.data(2,:)) - min(tsd.data(2,:)))/10),['Time: ' num2str(tsd.tvec(iF) - tsd.tvec(1),'%.2f%') '/' num2str(tsd.tvec(end) - tsd.tvec(1),'%.2f%')], 'color', 'w');
        %     drawnow
        pause(mode(diff(tsd.tvec))/100)
        
        %     clear p1
        delete(t1)
        delete(p1)
        delete(p2)
        
        
    end
    
end

%% outputs

boxes = [];
boxes.close = close_p;
boxes.open = open_p;
boxes.open_wR = open_wR_p;
boxes.open_wL = open_wL_p;

%% info

this_dir = cd;
parts = strsplit(this_dir, filesep); 

m_idx = 1+find(contains(parts, 'EPM')); 
d_idx = m_idx+2; 

info.subject = parts{m_idx};
info.date = parts{d_idx};



%% save the output if there is a save dir.
out.emp_idx = emp_idx;
out.labels = labels;
out.boxes = boxes;
out.dip_IV = dip_IV_m; 
out.CtT_IV = CtT_IV_m; 
out.OtT_IV = OtT_IV_m; 

if ~isempty(save_dir)
    save([save_dir filesep info.subject '_' info.date '.mat'], 'out')
end




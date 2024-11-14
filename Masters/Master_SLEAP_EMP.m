function Master_SLEAP_EMP(data_dir)

cd(data_dir) 

a_list = dir('*.avi');
fnames = dir('label*analysis.h5');

if length(a_list) ~= length(fnames)
    
    error('# of .avi files and # of label h5 files are not equal')
end



v_num =[]; 
for ii = 1:length(fnames)

    p_idx = strfind(fnames(ii).name, '.')+1;
    n_idx(1) = strfind(fnames(ii).name(p_idx(1):p_idx(end)), '_')+p_idx(1);
    n_idx(2) = strfind(fnames(ii).name, '.analysis.h5');
    v_num(ii) = str2double(fnames(ii).name(n_idx(1):n_idx(2)-1));
    
    lab_list{ii} = fnames(ii).name;
end

[~, s_idx] = sort(v_num); 

lab_list_s = lab_list(s_idx);

%% Load and scale the position data from Sleap

conv_f = [1 1]; 

tsd = MS_SLEAP2TSD(lab_list_s, 30, [],conv_f, 1);

if (tsd.tvec(end) - tsd.tvec(1)) > 600
   
    s_idx = nearest_idx(tsd.tvec(end) - 600, tsd.tvec);
    
    tsd.tvec = tsd.tvec(s_idx:end); 
    tsd.data = tsd.data(:,s_idx:end); 

    
end

%% plot to checl


plot(tsd.data(3,:), tsd.data(4,:), '.');


%% set the boarders and collect the 

imagesc(tsd.mean_frame)
closed = MS_drawrectangle_wait('Closed',[ 0 0 1]);
open_trk = MS_drawrectangle_wait('Open track',[1 1 0]);
open_wide_R = MS_drawrectangle_wait('Open wide right',[0 1 1]);
open_wide_L = MS_drawrectangle_wait('Open wide left',[1 0 1]);

[~, ~, close_p] = MS_rec2corner(closed.Position); 
[~, ~, open_p] = MS_rec2corner(open_trk.Position); 
[~, ~, open_wR_p] = MS_rec2corner(open_wide_R.Position); 
[~, ~, open_wL_p] = MS_rec2corner(open_wide_L.Position); 


xcross = intersect(open_p,close_p);

hold on
plot(xcross, 'FaceColor', [0 1 0], 'FaceAlpha', .5, 'EdgeColor', [0 1 0], 'LineWidth', 5); 

%% vectorize occupancy

head_idx = [3 4]; 
body_idx = [5 6]; %Which body node to use. 

figure(878)
clf
imagesc(tsd.mean_frame)
colormap('bone')

hold on

plot(close_p, 'FaceColor',[ 0 0 1], 'EdgeColor',[0 0 1], 'FaceAlpha', .2)
plot(open_p, 'FaceColor',[ 1 1 0], 'EdgeColor',[1 1 0], 'FaceAlpha', .2)
plot(xcross, 'FaceColor', [0 1 0], 'FaceAlpha', .5, 'EdgeColor', [0 1 0], 'LineWidth', 3); 
plot(open_wR_p, 'FaceColor', [0 1 1], 'FaceAlpha', .5, 'EdgeColor', [0 1 1], 'LineWidth', 3); 
plot(open_wL_p, 'FaceColor', [1 0 1], 'FaceAlpha', .5, 'EdgeColor', [1 0 1], 'LineWidth', 3); 


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

c_idx = c_idx & (~ x_idx & ~ o_idx);
o_idx = o_idx & (~ x_idx & ~ c_idx);
n_idx = ~x_idx & ~o_idx & ~c_idx;

% get putative head dips. 
[h_o_in, h_o_on] = inpolygon(tsd.data(head_idx(1),:), tsd.data(head_idx(2),:), open_p.Vertices(:,1), open_p.Vertices(:,2));

d_idx = ~(h_o_in | h_o_on)  & (~c_idx & ~ x_idx);


p1 = plot(tsd.data(body_idx(1),c_idx),tsd.data(body_idx(2),c_idx),'b+'); % points inside
p2 = plot(tsd.data(body_idx(1),o_idx),tsd.data(body_idx(2),o_idx),'yo'); % points outside
p3 = plot(tsd.data(body_idx(1),x_idx),tsd.data(body_idx(2),x_idx),'gs'); % points outside
p4 = plot(tsd.data(body_idx(1),n_idx),tsd.data(body_idx(2),n_idx),'kd'); % points outside
p5 = plot(tsd.data(head_idx(1),d_idx),tsd.data(head_idx(2),d_idx),'rx'); % points outside





lg = legend([p1 p2 p3 p4 p5], {'closed idx', 'open idx', 'x idx', 'none idx', 'dip idx'});
% lg.String(1:5)  = []; 
%% collect the data

emp_idx = NaN(size(tsd.data(3,:)));

emp_idx(n_idx) = 0;
emp_idx(c_idx) = 1;
emp_idx(o_idx) = 2;
emp_idx(x_idx) = 3;


fprintf('EMP: %0.1f%% open | %0.1f%% closed | %0.1f%% xcross | %0.1f%% none\n', (sum(o_idx)/length(o_idx))*100, (sum(c_idx)/length(c_idx))*100, (sum(x_idx)/length(x_idx))*100, (sum(n_idx)/length(n_idx))*100)

figure(8818)
clf

plot(emp_idx)
% set(gca,    % fill in states on y axis

% transitions
min_t = 2; % time in seconds for the animal to occuppy a region for it to count. 
min_t = min_t * round(1/mode(diff(tsd.tvec))); 



[c_entries] = find(diff(emp_idx == 1) > 0); 
[o_entries] = find(diff(emp_idx == 2) > 0); 
[x_entries] = find(diff(emp_idx == 3) > 0); 

[c_exits] = find(diff(emp_idx == 1) < 0); 
[o_exits] = find(diff(emp_idx == 2) < 0); 
[x_exits] = find(diff(emp_idx == 3) < 0); 

d_entries = find(diff(d_idx == 1) >0); 
d_exts    = find(diff(d_idx == 1) <0); 

% make 'trials' for each entry/exit type. 

if length(d_entries) > length(d_exits)
    d_exits(end+1) = length(d_idx); 
end



% use IV functions to link transitions or something. 







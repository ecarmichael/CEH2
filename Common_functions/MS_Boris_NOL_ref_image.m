function MS_Boris_NOL_ref_image(data_dir, env_diameter)


%% init

if nargin < 1
    
    data_dir = cd;
    env_diameter =  45; % in cm
elseif nargin < 2
    env_diameter = 45; 
end




%% grab the first frame from a video object


cd(data_dir)


% get the video

% v_name = dir('behav.avi'); 


v_obj = VideoReader('behav.avi')

f_1 = read(v_obj, 1); 


%% plot first frame
figure(6000)
clf

% imagesc()
a = imshow(f_1);


% try to auto find the objects

% [centersBright, radiiBright] = imfindcircles(f_1,[10 30], 'Sensitivity', .90,'ObjectPolarity','dark');

% viscircles(centersBright, radiiBright,'Color','b')
%% get the pixel to cm ratio, 



% ask for a line across the environment
line_obj = MS_drawline_wait;

box_dia = sqrt(((line_obj.Position(1,2)-line_obj.Position(1,1))^2)+((line_obj.Position(2,2)-line_obj.Position(2,1))^2));


px_ratio = box_dia / env_diameter; 

%% object 1
c_ord = MS_linspecer(3); 

obj_a = MS_drawcircle_wait;

obj_a.Visible = 'off'; 

viscircles(obj_a.Center, obj_a.Radius,'Color',c_ord(1,:))
% add in a second circle with an extra 2cm 

viscircles(obj_a.Center, obj_a.Radius+ (px_ratio*2),'Color',c_ord(1,:))

% REPEAT FOR obj 2


obj_b = MS_drawcircle_wait; 
obj_b.Visible = 'off'; 
viscircles(obj_b.Center, obj_b.Radius,'Color',c_ord(3,:))
% add in a second circle with an extra 2cm 

viscircles(obj_b.Center, obj_b.Radius+ (px_ratio*2),'Color',c_ord(3,:))

x_lim = xlim; 
y_lim = ylim; 

rectangle('Position',[x_lim(1), y_lim(1), x_lim(2) - x_lim(1), y_lim(2) - y_lim(1)],'edgecolor', 'k');
%% remove the og and save as a reference png


delete(a)


export_fig( 'NOL_ref.png', '-transparent', gcf)

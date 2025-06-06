function pattern_out = MS_radial_snap_grid(pos, bin_s, pattern)
%% MS_radial_snap_grid: defines a grid for the radial maze and snaps the arms to lines. If a 'pattern' is given as an input it will bypass the user feedback phase. 
%
%
%
%    Inputs: 
%    - pos: [struct]          positions data in the TSD format. should be in cm
%
%    - bin_s: [2x1 double]    size of x and y bins. Default is 3cm
%
%    - pattern: [struct]      output from MS_radial_snap_grid. If provided
%                              will data to this patter. 
%
%
%    Outputs: 
%    - pattern_out: [struct]      contains coordinates for isolating and
%    snapping position data to the radial grid. 
%
%
%
% EC 2023-08-30   initial version 
%
%
%
%% initialize


if nargin < 2
    bin_s = [ 3 3];
    pattern = [];
elseif nargin < 3
    pattern = [];
end

c_ord = linspecer(8); 


%% swap data from the miniscope format to tsd; 

pos_tsd = []; 
pos_tsd.data = pos.position' ;
pos_tsd.tvec = pos.time; 

pos = pos_tsd; 
%% create a basic grid
xmin = 20; % make sure these are bigger and consistent across sessions. 
xmax = 120; 
ymin = 0;
ymax = 80; 

x_edges = xmin:bin_s(1):xmax;
y_edges = ymin:bin_s(2):ymax;

occ_hist = histcn(pos.data(1:2,:)',x_edges,y_edges); % 2-D version of histc()


%% generate a pattern to snap to. 

figure(1)
clf
hold on
plot(pos.data(1,:), pos.data(2,:), '.')
MS_rescale_axis(pos.data(1,:), pos.data(2,:), 0.1)

axis equal

if isempty(pattern)

    cen_circle = drawcircle(gca,  'color', 'k'); 
    wait(cen_circle)
    
    pattern.c_cent = cen_circle.Center;
    pattern.c_rad = cen_circle.Radius; 
end

% convert circle to x and Y
t = 0:0.01:2*pi; % use to make a circle. 

c_x = cos(t)*pattern.c_rad+pattern.c_cent (1);
c_y = sin(t)*pattern.c_rad+pattern.c_cent (2);

plot(c_x, c_y, '.r')


c_idx = inpolygon(pos.data(1,:), pos.data(2,:), c_x, c_y);

pattern.c_idx = c_idx; 
plot(pos.data(1,c_idx), pos.data(2,c_idx), '.r')

%% draw rectangle and then rotate. 

this_axis = MS_drawrectangle_wait

%%
this_poly = polyshape([this_axis.Position(1), this_axis.Position(1) this_axis.Position(1)+this_axis.Position(3), this_axis.Position(1)+this_axis.Position(3)],...
    [this_axis.Position(2)+this_axis.Position(4), this_axis.Position(2),this_axis.Position(2), this_axis.Position(2)+this_axis.Position(4)]);
% plot(this_poly)

% NE_axis = this_axis; 
% NE_axis.RotationAngle = 45; 

% NE_poly = polyshape([NE_axis.Position(1), NE_axis.Position(1) NE_axis.Position(1)+NE_axis.Position(3), NE_axis.Position(1)+NE_axis.Position(3)],...
    % [NE_axis.Position(2)+NE_axis.Position(4), NE_axis.Position(2),NE_axis.Position(2), NE_axis.Position(2)+NE_axis.Position(4)]);

    NE_poly = rotate(this_poly, 45, pattern.c_cent);
plot(NE_poly)

% EW_axis = this_axis; 
% EW_axis.RotationAngle = 90; 
% 
% EW_poly = polyshape([EW_axis.Position(1), EW_axis.Position(1) EW_axis.Position(1)+EW_axis.Position(3), EW_axis.Position(1)+EW_axis.Position(3)],...
%     [EW_axis.Position(2)+EW_axis.Position(4), EW_axis.Position(2),EW_axis.Position(2), EW_axis.Position(2)+EW_axis.Position(4)]);
% 

    EW_poly = rotate(this_poly, 90, pattern.c_cent);

plot(EW_poly)

% SE_axis = this_axis; 
% SE_axis.RotationAngle = 135; 
% 
% SE_poly = polyshape([SE_axis.Position(1), SE_axis.Position(1) SE_axis.Position(1)+SE_axis.Position(3), SE_axis.Position(1)+SE_axis.Position(3)],...
%     [SE_axis.Position(2)+SE_axis.Position(4), SE_axis.Position(2),SE_axis.Position(2), SE_axis.Position(2)+SE_axis.Position(4)]);

    SE_poly = rotate(this_poly, 135, pattern.c_cent);

plot(SE_poly)

rad_poly = union(this_poly, SE_poly); 
rad_poly = union(rad_poly, EW_poly); 
rad_poly = union(rad_poly, NE_poly); 

%% plot again
hold on

plot(rad_poly, 'FaceColor','r')

%% try it with a line. 

    % north arm
if ~isfield(pattern, 'N')
    disp('Draw north arm')
    this_arm = drawline(gca, 'color', c_ord(1,:));
    pattern.N.pos = this_arm.Position; 

end
    

if ~isfield(pattern, 'NE')
        disp('Draw north east arm')
    this_arm = drawline(gca, 'color', c_ord(2,:));
    pattern.NE.pos = this_arm.Position; 

end

if ~isfield(pattern, 'E')
        disp('Draw east arm')
    this_arm = drawline(gca, 'color', c_ord(3,:));
    pattern.E.pos = this_arm.Position; 
end

% South East arm
if ~isfield(pattern, 'SE')
        disp('Draw south east arm')
   this_arm = drawline(gca, 'color', c_ord(4,:));
    pattern.SE.pos = this_arm.Position; 
end

% South  arm
if ~isfield(pattern, 'S')
        disp('Draw south arm')
   this_arm = drawline(gca, 'color', c_ord(5,:));
    pattern.S.pos = this_arm.Position; 
end
    

% South West arm
if ~isfield(pattern, 'SW')
        disp('Draw south west arm')
   this_arm = drawline(gca, 'color', c_ord(6,:));
   pattern.SW.pos = this_arm.Position; 
end
    

%  West arm
if ~isfield(pattern, 'W')
        disp('Draw west arm')
   this_arm = drawline(gca, 'color', c_ord(7,:));
    pattern.W.pos = this_arm.Position; 
end

%  North West arm
if ~isfield(pattern, 'NW')
        disp('Draw north west arm')
   this_arm = drawline(gca, 'color', c_ord(8,:));
    pattern.NW.pos = this_arm.Position; 
end

%% draw lines using the center of the maze and drawing a line out to the end from the center. 

figure(10)
clf
hold on
plot(pos.data(1,:), pos.data(2,:), '.')
MS_rescale_axis(pos.data(1,:), pos.data(2,:), 0.1)

axis equal

arms = {'N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'}; 

t = 0:0.01:2*pi; % use to make a circle. 

c_x = cos(t)*pattern.c_rad+pattern.c_cent(1);
c_y = sin(t)*pattern.c_rad+pattern.c_cent(2);

hold on
plot(c_x, c_y, '.r')
plot(pattern.c_cent(1), pattern.c_cent(2), 'x')

for iA = 1:length(arms)
    line([pattern.c_cent(1) pattern.(arms{iA}).pos(1,1)]', [pattern.c_cent(2) pattern.(arms{iA}).pos(1,2)]','color',  c_ord(iA,:), 'linewidth', 3)
    % line_grid(:, )
end




%%  Snap the arms using a polygon for each arm;
arms = {'N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'}; 

for iA = 1:length(arms)
    
    if isfield(pattern, arms{iA})
        
        if ~isfield(pattern.(arms{iA}), 'R_pos')
            
            fprintf('Draw rectangle for <strong>%s</strong> arm. put the start of the rectangle at the end. Rotate as needed\n', arms{iA})
            this_arm = drawpolygon(gca, 'color', c_ord(iA,:));
            fprintf('Press enter when done to move to next arm\n')

            wait(this_arm)
            
            pattern.(arms{iA}).R_pos = this_arm.Position;
%             pattern.(arms{iA}).R_ang = this_arm.RotationAngle;


% %% draw a poly based on the line. 
% 
%             x = pattern.(arms{iA}).pos(:,1);
%             y = pattern.(arms{iA}).pos(:,2);
%             
%             h_x = x;  
%             h_y = [y(1) y(1)]; 
%             
%             adj = x(2) - x(1); 
%             opp = y(2) - y(1);
%             
% %             theta = tan(opp/adj); 
%             
% %             Theta = atan2(norm(cross((, v_2)), dot(v_1, v_2));
% %             r_pos_x = [
%             
            
        end
    end
end

%% get inclusion idx for each arm. 

for iA = 1:length(arms)
    
    % get the indicies within the arm rectangle. 
    this_r = pattern.(arms{iA}).R_pos; 
    [in_idx, on_idx]= inpolygon(pos.data(1,:), pos.data(2,:), this_r(:,1), this_r(:,2)); 
    %     [in_idx, on_idx]= inpolygon(pos.data(1,:), pos.data(2,:), [this_r(1), this_r(1), this_r(1)+this_r(3), this_r(1)+this_r(3)], [this_r(2),this_r(2)+this_r(4),this_r(2)+this_r(4), this_r(2)] );

    
    pattern.(arms{iA}).idx = (in_idx | on_idx) & ~c_idx; 
    
    
    plot(pos.data(1, pattern.(arms{iA}).idx),  pos.data(2, pattern.(arms{iA}).idx), '.', 'color', c_ord(iA,:))
    
end


%% Check for points that fall outside of the center or arms. 

inc_idx = pattern.c_idx | pattern.N.idx | pattern.NE.idx | pattern.E.idx | pattern.SE.idx | pattern.S.idx | pattern.SW.idx | pattern.W.idx | pattern.NW.idx; 

% fill in missing data. 

pos.data(1,~inc_idx) = NaN; 
pos.data(2,~inc_idx) = NaN; 


pos.data(1,:) = fillmissing(pos.data(1,:), 'nearest');
pos.data(2,:) = fillmissing(pos.data(2,:), 'nearest');

%% Snap to line in arm. 
linpos_temp(keep_idx) = griddata(Coord_in.coord(1,:),Coord_in.coord(2,:),coord_vals,x(keep_idx),y(keep_idx),'nearest');

%% try to make a hex grid
% 
% Rad3Over2 = sqrt(3) / 2;
% [X Y] = meshgrid(0:2:xmax+3);
% n = size(X,1);
% X = Rad3Over2 * X;
% Y = Y + repmat([0 0.5],[n,n/2]);
% 
% Plot the hexagonal mesh, including cell borders
% [XV YV] = voronoi(X(:),Y(:)); plot(XV,YV,'b-')
% axis equal, axis([10 20 10 20]), zoom on
% 
% 
% histogram2(pos.data(1,:), pos.data(2,:), XV(1,:), YV(2,:))
end



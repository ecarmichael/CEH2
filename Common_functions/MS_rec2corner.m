function [x,y, poly_out] = MS_rec2corner(position); 

% disp(position)


x = [position(1), position(1), position(1)+position(3), position(1)+position(3)];
y = [position(2), position(2)+position(4),  position(2)+position(4),position(2)];

poly_out = polyshape(x, y); 
function MS_imagesc_append_values(M, x, y, colour)

[rows, cols] = size(M);

if nargin < 2
    x = 1:cols; 
    y = 1:rows;
    colour = [1 1 1]; 
elseif  nargin < 3
    y = 1:rows; 
    colour = [1 1 1];
elseif nargin < 4
    colour = [1 1 1]; 
end


for r = 1:rows
    for c = 1:cols
        % Extract value and convert to string
        valStr = num2str(M(r, c), 3); 
        
        % Place text (Note: text function expects X=column, Y=row)
        text(x(c), y(r), valStr, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'Color', colour, ...
            'FontSize', 12, ...
            'FontWeight', 'bold');
    end
end
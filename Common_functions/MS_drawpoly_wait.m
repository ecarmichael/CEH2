function hROI =  MS_drawpoly_wait(label, c_ord)
%% MS_drawrectangle_wait:  draw a rectangle and wait for the user to finish before gettting the ROI.
% from: https://www.mathworks.com/help/images/use-wait-function-after-drawing-roi-example.html

fprintf('<strong>Draw a polygone and then double click on the shape to pass the function and get the roi</strong>\n')

if nargin < 1
    label = ' ';
    c_ord = [0 0 0]; 
end

hROI = drawpolygon('Label',label,'Color',c_ord); 

l = addlistener(hROI,'ROIClicked',@clickCallback);

% Block program execution
uiwait;

% Remove listener
delete(l);

end

function clickCallback(~,evt)

if strcmp(evt.SelectionType,'double')
    uiresume;
end

end
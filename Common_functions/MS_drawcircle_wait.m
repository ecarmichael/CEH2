function hROI =  MS_drawcircle_wait
%% MS_drawcircle_wait:  draw a cirle and wait for the user to finish before gettting the ROI.

% from: https://www.mathworks.com/help/images/use-wait-function-after-drawing-roi-example.html

fprintf('<strong>Draw a circle and then double click on the shape to pass the function and get the roi</strong>\n')

hROI = drawcircle; 

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
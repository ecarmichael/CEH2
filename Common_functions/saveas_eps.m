function saveas_eps(fname, save_dir)
%% saveas_eps: this is a workaround to make sure that matlab saves an eps
%file that is trimmed to the correct size
%
% uses the method outlined here: https://www.mathworks.com/matlabcentral/answers/162283-why-is-the-figure-in-my-eps-file-generated-using-matlab-r2014b-in-the-wrong-position-and-with-extra
%
%  inputs: h : figure hangle(can be "gcf")
%          fname: [string] contains the name and the path to the save
%          directory
%

h = get(gcf);
D = h.PaperPosition; % Returns 1x4 vector [left right width height]
h.PaperSize = [D(3) D(4)]; %default PaperSize is [8.5 11]
    pushdir(save_dir);
    eval(sprintf('print -depsc2 -tiff  -r300 -painters %s',fname));
    popdir;

end
%%
close all

maps_e = MS_get_place_field([], restrict(out.Encode.S, out.Encode.trials(1,:), out.Encode.trials(2,:)), restrict(out.Encode.pos, out.Encode.trials(1,:), out.Encode.trials(2,:)), out.Encode.speed)

pause(1)
set(102, 'position', [ 2000 50 1000 800])
pause(1)
figure(102)
exportgraphics(gcf, ['/home/williamslab/Williams Lab Dropbox/Eric Carmichael/CIHR_2023_September' filesep 'Example_PCs_radialMaze' filesep out.meta.subject '_' out.meta.session  '_place_encode_1.pdf'], 'ContentType', 'vector');
close(102)

set(103, 'position', [ 2000 50 1000 800])
pause(1)
figure(103)
exportgraphics(gcf, ['/home/williamslab/Williams Lab Dropbox/Eric Carmichael/CIHR_2023_September' filesep 'Example_PCs_radialMaze' filesep out.meta.subject '_' out.meta.session  '_place_encode_2.pdf'], 'ContentType', 'vector');

%%
close all

maps_r = MS_get_place_field([], restrict(out.Recall.S, out.Recall.trials(1,:), out.Recall.trials(2,:)), restrict(out.Recall.pos, out.Recall.trials(1,:), out.Recall.trials(2,:)), out.Recall.speed)


pause(1)
set(102, 'position', [ 2000 50 1000 800])
pause(1)
figure(102)
exportgraphics(gcf, ['/home/williamslab/Williams Lab Dropbox/Eric Carmichael/CIHR_2023_September' filesep 'Example_PCs_radialMaze' filesep out.meta.subject '_' out.meta.session  '_place_recall_1.pdf'], 'ContentType', 'vector');
close(102)

pause(1)
set(103, 'position', [ 2000 50 1000 800])
pause(1)
figure(103)

exportgraphics(gcf, ['/home/williamslab/Williams Lab Dropbox/Eric Carmichael/CIHR_2023_September' filesep 'Example_PCs_radialMaze' filesep out.meta.subject '_' out.meta.session  '_place_recall_2.pdf'], 'ContentType', 'vector');
close(103)

close all
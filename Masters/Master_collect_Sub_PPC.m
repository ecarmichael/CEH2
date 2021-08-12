function Master_collect_Sub_PPC(inter_dir)

c_ord = linspecer(5);

cd(inter_dir)

Blocks = {'Pre_sleep', 'W_maze', 'OF', 'Post_sleep'};
%% look for all the PPC output files


fnames = dir('*.mat');

for ii =  1:length(fnames)
    if strcmp(fnames(ii).name(1:3), 'PPC')
        ppc_files{ii} = fnames(ii).name;
    end
end

ppc_files(cellfun('isempty', ppc_files)) = [];


%% check for incomplete files

keep_idx = zeros(1,length(ppc_files)); 
for iP = length(ppc_files):-1:1
    load(ppc_files{iP}, 'PPC');
    
    if ~isfield(PPC, Blocks{4}) || ~isfield(PPC.(Blocks{4}), 'z_ppc')
        disp(ppc_files{iP})
        keep_idx(iP) = 0;
    else
        keep_idx(iP) = 1;
    end
end

ppc_files(~keep_idx) = [];

%% loop over files and collect the z_ppc and file name. 

z_mat = [];
for iP = length(ppc_files):-1:1
    load(ppc_files{iP}, 'PPC');

    % collect the names and z_ppc
    z_mat.label{iP} = ppc_files{iP}(4:end-4);
    
    z_mat.cell{iP} = ppc_files{iP}(strfind(ppc_files{iP}, 'TT'): strfind(ppc_files{iP}, 'TT')+7);
    z_mat.subject{iP} = ppc_files{iP}(4:6);
    z_mat.date{iP} = ppc_files{iP}(8:strfind(ppc_files{iP}, 'TT')-2);
    z_mat.day{iP} = z_mat.date{iP}(end-1:end);
    for iB =1:length(Blocks)
       z_mat.(Blocks{iB})(iP,:) = PPC.(Blocks{iB}).z_ppc; 
    end
    z_mat.obs_freq = PPC.Pre_sleep.obs_freq; 
    
    clear PPC
    
end

%% make a few plots
close all
figure('PaperPositionMode', 'auto')
for iB =1:length(Blocks)
   subplot(1,4,iB)
   imagesc(z_mat.obs_freq, 1:length(z_mat.cell),  z_mat.(Blocks{iB}));
   xlabel('frequency (hz)')
   ylabel('cell #')
   if strcmp(Blocks{iB}, 'OF')
          title('Open Field');
   else
   title(strrep(Blocks{iB}, '_', ' '));
   end
   caxis([1.98 8])
   set(gca, 'ytick', 1:length(z_mat.cell));
   set(gca, 'Box', 'off')
   ylim([.25 length(z_mat.cell)])
   rectangle('position',[1, .25, 3, .25], 'facecolor', c_ord(5,:), 'edgecolor', [0 0 0 0])
   rectangle('position',[6, .25, 6, .25], 'facecolor', c_ord(2,:), 'edgecolor', [0 0 0 0])
%    rectangle('position',[6, .25, 6, .25], 'facecolor', c_ord(5,:), 'edgecolor', [0 0 0 0])
   rectangle('position',[40, .25, 40, .25], 'facecolor', c_ord(1,:), 'edgecolor', [0 0 0 0])
   rectangle('position',[80, .25, 21, .25], 'facecolor', c_ord(3,:), 'edgecolor', [0 0 0 0])

end
cb=colorbar;
cb.Position(1) = cb.Position(1) + .08; 
cb.Position(3) = cb.Position(3) * 2;
cb.Label.String = 'z-score';
colormap([0.2,0.2,0.2; parula(512)]);
SetFigure([], gcf)
set(gcf, 'position', [70 247 1748 638]);

print('PPC_summary', '-dpng', '-r300')
saveas(gcf, 'PPC_summary', 'fig')
print('PPC_summary', '-dsvg'); 

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters PPC_summary.pdf

% print(gcf, '-dpdf', 'PPC_summary.pdf');
% saveas(gcf, 'PPC_summary', 'svg')

% 
% %% same thing but vertical
% figure(102)
% for iB =1:length(Blocks)
%    subplot(4,1,iB)
%    imagesc(z_mat.obs_freq, 1:length(z_mat.cell),  z_mat.(Blocks{iB}));
%    ylabel('cell #')
%    title(strsplit(Blocks{iB}, '_'), 'position', [-10 median(1:length(z_mat.cell))]*1.3)
%    caxis([1.98 8])
%    if iB == length(Blocks)
%        xlabel('frequency (hz)')
%    else
%        set(gca, 'xtick', []);
%    end
%    set(gca, 'ytick', 1:length(z_mat.cell));
%    set(gca, 'Box', 'off')
%    ylim([.25 length(z_mat.cell)])
%    rectangle('position',[1, .25, 3, .25], 'facecolor', c_ord(5,:), 'edgecolor', [0 0 0 0])
%    rectangle('position',[6, .25, 6, .25], 'facecolor', c_ord(2,:), 'edgecolor', [0 0 0 0])
% %    rectangle('position',[6, .25, 6, .25], 'facecolor', c_ord(5,:), 'edgecolor', [0 0 0 0])
%    rectangle('position',[40, .25, 40, .25], 'facecolor', c_ord(1,:), 'edgecolor', [0 0 0 0])
%    rectangle('position',[80, .25, 21, .25], 'facecolor', c_ord(3,:), 'edgecolor', [0 0 0 0])
% end
% set(gcf, 'position', [600 50 784 588]);
% 
% cb=colorbar;
% cb.Position(1) = cb.Position(1) + .05; 
% cb.Position(3) = cb.Position(3) * 2;
% cb.Label.String = 'zscore';
% 
% SetFigure([], gcf)
% % 
% % saveas(gcf, 'PPC_summary', 'png')
% % saveas(gcf, 'PPC_summary', 'fig')
% % print('PPC_sumary', '-depsc');


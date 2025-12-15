%% MASTER_Asmbly_debug

restoredefaultpath
usr = char(java.lang.System.getProperty('user.name')); 

if strcmp(computer, 'GLNXA64')
    codebase_dir = ['/home/' usr '/Documents/Github/vandermeerlab/code-matlab/shared'];
    ca_dir = ['/home/' usr '/Documents/Github/CEH2'];
    % oasis_dir = ['/home/' usr '/Documents/Github/OASIS_matlab'];
    code_dir = ['/home/' usr '/Williams Lab Dropbox/Eric Carmichael/Comp_Can_inter/Git_repos/Dos-Santos Assembly ICA'];
    % RnR_dir = ['/home/' usr '/Williams Lab Dropbox/Eric Carmichael/Comp_Can_inter/Git_repos/RnR_methods'];
    main_dir = ['/home/' usr '/'];
    
elseif strcmp(computer, 'MACA64')
    codebase_dir = ['/Users/' usr '/Documents/Github/vandermeerlab/code-matlab/shared'];
    ca_dir = ['/Users/' usr '/Documents/Github/CEH2'];
    % oasis_dir = ['/Users/' usr '/Documents/Github/OASIS_matlab'];
    code_dir = ['/Users/' usr '/Williams Lab Dropbox/Eric Carmichael/Comp_Can_inter/Git_repos/Dos-Santos Assembly ICA'];
    % RnR_dir = ['/Users/' usr '/Williams Lab Dropbox/Eric Carmichael/Comp_Can_inter/Git_repos/RnR_methods'];
    main_dir = ['/Users/' usr filesep];
else
    codebase_dir = ['C:\Users\' usr '\Documents\Github\vandermeerlab\code-matlab\shared'];
    ca_dir = ['C:\Users\' usr '\Documents\Github\CEH2'];
    % oasis_dir = ['C:\Users\' usr '\Documents\Github\OASIS_matlab'];
    code_dir = ['C:\Users\' usr '\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\Git_repos\Dos-Santos Assembly ICA'];
    % RnR_dir = ['C:\Users\' usr '\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\Git_repos\RnR_methods'];
    main_dir = ['C:\Users\' usr '\'];
end

restoredefaultpath
c_d = cd;

addpath(genpath(ca_dir));
addpath(genpath(codebase_dir))
% addpath(genpath(RnR_dir));
addpath(genpath(code_dir))

cd(c_d)

move_thresh  = 9;
bin_size = .5;

inter_dir = strrep([main_dir 'Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\PC9'], '\', filesep);
fig_dir = [inter_dir filesep 'checks_proof'];

%%  Extract assembly data
cd(inter_dir);

f_list = dir('*data*');

B_out = cell(length(f_list), 1); 
session = B_out; 
novel_idx = []; anx_idx = []; HS_idx = [];
method = 'binary';


for ii = 1:length(f_list)
    session{ii} = f_list(ii).name;
    % compute assemblies and related ReActs
    % A_out{ii} = Pipeline_Asmbly(f_list(ii).name,bin_size, move_thresh, method);
    % P_out{ii} = Pipeline_Asmbly_place(f_list(ii).name,bin_size, move_thresh, method);

    B_out{ii} = Pipeline_Asmbly_top_cells(f_list(ii).name,bin_size, move_thresh, method, [], 256);

    B_out{ii} = Pipeline_Asmbly_append_SWS(f_list(ii).name, B_out{ii});

    B_out{ii} = Pipeline_Asmbly_append_preA(B_out{ii});

    % Summary plots
    %                 Pipline_Asmbly_plot(A_out{ii}, [fig_dir filesep method]);
    %     Pipline_Asmbly_plot(P_out{ii}, [fig_dir filesep method filesep 'place']);

    % Pipline_Asmbly_plot(B_out{ii}, [fig_dir filesep 'bin_' num2str(bin_size) filesep method filesep 'best']);

    % Pipline_Asmbly_plot_SWS(B_out{ii}, [fig_dir filesep 'bin_' num2str(bin_size) filesep method filesep 'best_SWS']);
    
    close all
    
    if ~isempty(strfind(f_list(ii).name, 'HATDS'))
        HS_idx(ii) = 1;
    else
        HS_idx(ii) = 0;
    end
    
    if ~isempty(strfind(f_list(ii).name, 'D1')) %|| ~isempty(strfind(f_list(ii).name, 'HATDS'))
        novel_idx(ii) = 1;
    end
    
    if ~isempty(strfind(f_list(ii).name, 'D5'))
        novel_idx(ii) = 0;
    end
    
    if ~isempty(strfind(f_list(ii).name, 'HAT'))
        anx_idx(ii) = 1;
    else
        anx_idx(ii) = 0;
    end
    
    if anx_idx(ii)== 0 && novel_idx(ii)==1
        %         B_out{ii} = Pipeline_Asmbly_append_preA(B_out{ii});
    end
    
end
novel_idx = logical(novel_idx);
anx_idx = logical(anx_idx);
HS_idx = logical(HS_idx);

%% save the assembly structure
save([inter_dir  strrep('\Assembly\inter\B_out_', '\', filesep) method '.mat'], 'B_out')


%% load data that has been processed.

method = 'binary';

load([inter_dir  strrep('\Assembly\inter\B_out_', '\', filesep) method '.mat'], 'B_out')

% rename due to progressive structure naming convention above.
A_out = B_out;
clearvars('B_out')

% remove pv1254 due to inconsistent running.
exclude_mouse = {'pv1254'};

for ii = length(A_out):-1:1
    if contains(A_out{ii}{1}.info.subject, exclude_mouse)
        fprintf('Removing sesson: <strong>%s</strong>\n', A_out{ii}{1}.info.subject);
        rm_idx(ii) = true;
    else
        rm_idx(ii) = false;
    end
end
A_out(rm_idx) = [];


novel_idx = []; anx_idx = []; HS_idx = [];
for iA = size(A_out,2):-1:1
    this_f = A_out{iA}{1}.info.session;
    
    %A_out{iA} = Pipeline_Asmbly_append_preA(A_out{iA});
    %A_out{iA} = Pipeline_Asmbly_append_postA(A_out{iA});
    
    if contains(this_f, 'HATDS')
        HS_idx(iA) = 1;
    else
        HS_idx(iA) = 0;
    end
    
    if contains(this_f, 'D1') %|| ~isempty(strfind(f_list(ii).name, 'HATDS'))
        novel_idx(iA) = 1;
    end
    
    if contains(this_f, 'D5')
        novel_idx(iA) = 0;
    end
    
    if contains(this_f, 'HAT')
        anx_idx(iA) = 1;
    else
        anx_idx(iA) = 0;
    end
    
end
function A_in = Pipeline_Asmbly_append_postA(A_in)
%% check the pre_rem for assemblies, remove assemblies without positive weights, apply assemblies to the wake track data

for ii = 1:length(A_in)
    if isfield(A_in{ii}, 'REM_temp'); A_in{ii} = rmfield(A_in{ii}, 'REM_temp'); end
    if isfield(A_in{ii}, 'REM_Wake_proj'); A_in{ii} = rmfield(A_in{ii}, 'REM_Wake_proj'); end
    if isfield(A_in{ii}, 'REM_A_pos'); A_in{ii} = rmfield(A_in{ii}, 'REM_A_pos'); end
    if isfield(A_in{ii}, 'pREM_Place_map'); A_in{ii} = rmfield(A_in{ii}, 'pREM_Place_map'); end
end

%% get the pre REM templates
for ii = 1:length(A_in)

    if size(A_in{ii}.REM_Post_data, 1) < size(A_in{ii}.REM_Post_data, 2)
        fprintf('<strong>Pre REM data has less samples (%0.0f) than cells (%0.0f)</strong> \n',size(A_in{ii}.REM_Pre_data, 1), size(A_in{ii}.REM_Pre_data, 2))
        A_in{ii}.postREM_temp = [];
        A_in{ii}.postREM_proj = [];
        A_in{ii}.postREM_Wake_proj = [];
        A_in{ii}.postREM_A_pos = [];
        continue
    else

        REM_temp_all = assembly_patterns(A_in{ii}.REM_Post_data);

        postREM_proj = assembly_activity(REM_temp_all, A_in{ii}.REM_Post_data');


        Wake_proj = assembly_activity(REM_temp_all, A_in{ii}.wake_data');

        [REM_temp, Wake_proj, REM_pos]=  MS_Asmbly_select(REM_temp_all, Wake_proj, 2);

        [~, postREM_proj, ~]=  MS_Asmbly_select(REM_temp_all, postREM_proj, 2);

        A_in{ii}.postREM_temp = REM_temp;
        A_in{ii}.postREM_proj = postREM_proj;
        A_in{ii}.postREM_Wake_proj = Wake_proj;
        A_in{ii}.postREM_A_pos = REM_pos;

        fprintf('[%.0f/%.0f = %.0f%%] Post REM Assemblies had cells with positive weights (%.2fs binsize)\n',size(REM_temp,2),size(REM_temp_all,2),  (size(REM_temp,2)/size(REM_temp_all,2))*100, A_in{ii}.info.bin)

        [A_in{ii}.postREM_stats, A_in{ii}.postREM_shuff.data, A_in{ii}.postREM_shuff.proj] = MS_Asmbly_proj_thresh(A_in{ii}.REM_Post_data, REM_temp, 500, 99);

        % rng(123,'twister')
        % nShuff = 100;
        % wake_shuff_mat = [];
        %
        % Ass_shuff = NaN(1,nShuff);
        % for iS = nShuff:-1:1
        %     tic
        %     shuff_data = NaN(size(A_in{ii}.REM_Pre_data));
        %     for ic = 1:size(A_in{ii}.REM_Pre_data,2)
        %         shuff_data(:,ic) = circshift(A_in{ii}.REM_Pre_data(:,ic), floor(MS_randn_range(1,1,1,size(A_in{ii}.REM_Pre_data,1))));
        %     end
        %
        %     this_ass = assembly_patterns(shuff_data');
        %     if ~isempty(this_ass)
        %         S_prog = assembly_activity(this_ass,shuff_data');
        %
        %         wake_shuff_mat(iS,:) =  S_prog(1,:);
        %         keep_idx(iS) = 1;
        %     else
        %         wake_shuff_mat(iS,:) = NaN;
        %         keep_idx(iS) = 0;
        %     end
        %     %     for ii = size(this_ass,2):-1:1
        %
        %     if sum(max(this_ass) > 0.2) >0
        %         Ass_shuff(iS) = sum(max(this_ass) > 0.2);
        %     else
        %         Ass_shuff(iS) = 0;
        %     end
        %     %     end
        %     fprintf('Shuff # %.0f found %.0f assemblies and took %2.2f seconds\n', iS, size(this_ass,2), toc)
        % end
        %
        % W_threshold = prctile(wake_shuff_mat(wake_shuff_mat >0), 99, 'all');
    end
end

%% get the mean place field templates;

min_N_place = 3;
for iB = length(A_in):-1:1
    if isempty(A_in{ii}.postREM_A_pos)
        A_in{iB}.postREM_Place_map = [];
        continue
    else
        %     if
        [map_out, place_idx] = MS_Asmbly_map(A_in{iB}.postREM_A_pos, A_in{iB}.place, min_N_place);

        Place_map = map_out;
        Place_map(~place_idx) = [];

        A_in{iB}.postREM_Place_map = Place_map;
    end
    %     Place_temp{iB} = P_temp{iB}(:,place_idx{iB});
    %     Place_proj{iB} = P_proj{iB}(:,place_idx{iB});
    %
    %     fprintf('[%.0f/%.0f = %.0f%%] Assemblies contained at least %0.0f place cells (%.2fs binsize)\n',size(Place_temp{iB},2),size(A_temp{iB},2),  (size(Place_temp{iB},2)/size(A_temp{iB},2))*100, min_N_place, bin_s(iB))
    %
end
%%
win_s = 2;
thresh = 10;
for iB = length(A_in):-1:1
    [A_in{ii}.postREM_wake_P_loc] = MS_Asmbly_act_loc(A_in{ii}.postREM_Wake_proj, A_in{ii}.wake_tvec, A_in{ii}.behav, win_s, thresh, 2/A_in{1}.bins);
end

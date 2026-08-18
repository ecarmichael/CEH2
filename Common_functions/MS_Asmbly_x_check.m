function [ReAct_out] = MS_Asmbly_x_check(A_out, ref, target)


rng(123, 'twister'); % for reproducibility.
opts = [];
opts.threshold.method = 'MarcenkoPastur';
opts.Patterns.method = 'ICA';
opts.Patterns.number_of_iterations = 500;
opts.threshold.permutations_percentile= 95;
opts.threshold.number_of_permutations= 500;

% init vars for number of assemblies per session and the shuff stats
ReAct_out.J_n_ass = [];  ReAct_out.J_r_ass = [];

for iA = length(A_out):-1:1

    % current assembly wieghts;
    if strcmpi(ref, 'pre')
        ref_temp = A_out{iA}{1}.pREM_temp;
    elseif strcmpi(ref, 'wake')
        ref_temp = A_out{iA}{1}.P_temp;
    elseif strcmpi(ref, 'post')
        ref_temp = A_out{iA}{1}.postREM_temp;
    end

    % target
    if strcmpi(target, 'pre')
        ref_tvec = A_out{iA}{1}.REM_Pre_tvec;
    elseif strcmpi(target, 'wake')
        ref_tvec = A_out{iA}{1}.wake_tvec; 
    elseif strcmpi(target, 'post')
        ref_tvec = A_out{iA}{1}.REM_Post_tvec; 
    end

        ReAct_out.sess_id{iA}  = [A_out{iA}{1}.info.subject '-' A_out{iA}{1}.info.session];

if isempty(ref_temp)
    ReAct_out.J_n_ass(iA, :) = NaN(1,length(A_out));
    ReAct_out.J_r_ass(iA, :) = NaN(1,length(A_out));
    ReAct_out.S_n_ass(iA, :) = NaN(1,length(A_out));
    ReAct_out.S_r_ass(iA, :) = NaN(1,length(A_out));

    continue
end

    % loop over sessions
    for jj = length(A_out):-1:1

        % test if A_temps can be found in J-data
        if strcmpi(target, 'pre')
            data_name = 'REM_Pre_data';
        elseif strcmpi(target, 'wake')
            data_name = 'wake_data'; 
        elseif strcmpi(target, 'post')
            data_name = 'REM_Post_data'; 
        end

        fprintf('Ref: %s - Target %s  ', ref, data_name)

        rng(123, 'twister'); % for reproducibility.
        A_alt_proj = assembly_activity(ref_temp,A_out{jj}{1}.(data_name)');

        rng(123, 'twister'); % for reproducibility.
        A_ref_proj = assembly_activity(ref_temp,A_out{iA}{1}.(data_name)');



        % trim the J_proj to be the same length as the A_proj
        if length(A_alt_proj) > length(A_ref_proj)
            A_alt_proj = A_alt_proj(:,1:length(A_ref_proj));
        end

        % alternative session
        rng(123, 'twister'); % for reproducibility.
        [Ref_stats, shuff.data, shuff.proj] = MS_Asmbly_proj_thresh(A_out{jj}{1}.(data_name), ref_temp, 500, 99);
        Ref_stats.p_val = [];
        Ref_stats.rate = [];
        Ref_stats.rate_p = [];
        Ref_stats.shuff_n = [];
        Ref_stats.shuff_r = [];

        for ii = size(A_alt_proj,1):-1:1
            Ref_stats.p_val(ii) = sum(sum(shuff.data > Ref_stats.R_thresh,2) > sum(A_alt_proj(ii,:) > Ref_stats.R_thresh))/ size(shuff.data,1);
            Ref_stats.rate(ii) = sum(A_alt_proj(ii,:) > Ref_stats.R_thresh) / ((ref_tvec(end) - ref_tvec(1))/60);
            Ref_stats.shuff_rate = sum(shuff.data > Ref_stats.R_thresh,2)./ ((ref_tvec(end) - ref_tvec(1))/60);
            Ref_stats.rate_p(ii) = sum(Ref_stats.shuff_rate > Ref_stats.rate(ii)) / length(Ref_stats.shuff_rate);

            Ref_stats.shuff_n(ii) = mean(sum(shuff.data > Ref_stats.R_thresh,2));
            Ref_stats.shuff_r(ii) = mean(sum(shuff.data > Ref_stats.R_thresh,2) / ((ref_tvec(end) - ref_tvec(1))/60));

        end

        ReAct_out.J_n_ass(iA, jj) = sum(Ref_stats.p_val < 0.05); % number of jj assemblies passing the pval test
        ReAct_out.J_r_ass(iA, jj) = mean(Ref_stats.rate(Ref_stats.rate_p < 0.05)); % number of jj assemblies passing the RATE pval test
        ReAct_out.S_n_ass(iA, jj) = mean(Ref_stats.shuff_n); % 
        ReAct_out.S_r_ass(iA, jj) = mean(Ref_stats.shuff_r); % 

        fprintf('iA: %d  | jj: %d\n', iA, jj)

    end

end


disp('done')
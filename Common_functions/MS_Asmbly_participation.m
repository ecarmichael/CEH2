function [p, p_in_mem, p_in_n_mem] = MS_Asmbly_participation(Proj_in, Temp_in, data_in, R_thresh,method, plot_flag)
%% MS_Asmbly_participation:
%
%
%
%    Inputs:
%    -
%
%
%
%    Outputs:
%    -
%
%
%
%
%  2024-07-02   initial version
%
%
%
%% initialize

if nargin < 5
    method = 'SD';
    plot_flag = false;
elseif nargin < 6
    plot_flag = false;
end

fprintf('Finding members using: <strong>%s</strong>\n', method)


p = cell(1, size(Temp_in, 2));
p_in_mem = nan(1, size(Temp_in, 2));
p_in_n_mem = nan(1, size(Temp_in, 2));


p = [];
%% loop over assemblies and get the

for iA = size(Temp_in,2):-1:1

    % grab the 'significant cells in the assembly
    z_weight = zscore(Temp_in(:,iA));

    if strcmpi(method, 'otsu')

        [o_data] = abs(Temp_in(:,iA));
        o_thresh = otsuthresh(o_data);
        o_thresh = o_thresh*(max(o_data)-min(o_data)) + min(o_data);
        M_idx = Temp_in(:,iA) > o_thresh;
    else

        M_idx = z_weight > 1.96; % pull out the indicies of the best cells.
    end
    % identify significant reactivations.
    R_idx = find(Proj_in(iA,:) > R_thresh);

    % check the data for the cell by cell activity across each reactivation
    A_part= [];

    if isempty(R_idx)  %|| length(R_idx) < 3
        continue
    else
        for iR = length(R_idx):-1:1

            this_R = data_in(R_idx(iR),:)>0';
            this_R(this_R & M_idx') =2;
            A_part(iR,:) = this_R';
        end


        % get the probability of being active in the assembly vs outside
         p{iA}.mem_idx = M_idx; 
        p{iA}.in_A = sum(A_part, 1) ./ size(A_part, 1);

        p{iA}.in_A_mem = mean(p{iA}.in_A(M_idx));
        p{iA}.in_A_n_mem = mean(p{iA}.in_A(~M_idx));

        p{iA}.out_A = sum(data_in(Proj_in(iA,:) < R_thresh,:)>0',1) ./ sum(Proj_in(iA,:) < R_thresh);

        % keep single value per assemblies for simple plotting/sorting
        % Store the computed probabilities in the output variables
        p_in_mem(iA) = p{iA}.in_A_mem;
        p_in_n_mem(iA) = p{iA}.in_A_n_mem;
        
        % Optionally, display the results in the command window
        if plot_flag
            disp(['Assembly ', num2str(iA), ': Mem = ', num2str(p{iA}.in_A_mem), ', Non-Mem = ', num2str(p{iA}.in_A_n_mem)]);
        end
        if plot_flag
            figure(iA); clf
            subplot(1,8,1)
            imagesc(z_weight);

            subplot(1,8,2:5)
            imagesc(A_part')

            subplot(1,8,6:7)
            imagesc(1:2,1:length(p{iA}.in_A),  [p{iA}.in_A; p{iA}.out_A; p{iA}.in_A.*p{iA}.out_A]');
        end
    end




end




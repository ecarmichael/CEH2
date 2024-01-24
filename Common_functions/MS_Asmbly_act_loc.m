function   [out] = MS_Asmbly_act_loc(P_proj,wake_tvec, behav, win_in, thresh, min_dist)
%% MS_Asmbly_act_loc:
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
% EC 2024-01-24   initial version
%
%
%
%% initialize
out = [];
win = floor(win_in * mode(diff(behav.time)));


%%
for ii = size(P_proj,1):-1:1
    
    [p_val, p_idx] = findpeaks(P_proj(ii,:),'MinPeakHeight', thresh ,'MinPeakDistance', min_dist);
    
    out{ii} = [];
    for ip = 1:length(p_idx)
        this_idx = nearest_idx(wake_tvec(p_idx(ip)), behav.time/1000);
        
        out{ii}.loc_time(ip,:) = behav.time(this_idx);
        
        out{ii}.loc(ip) = behav.position(this_idx,1)';
        
        if ((this_idx - win) >=1) && ((this_idx + win)<= length(behav.time))
            out{ii}.loc_mat(ip,:) = behav.position(this_idx - win:this_idx+win,1);
            
        else
            
            out{ii}.loc_mat(ip,:) = NaN((win*2)+1,1);
        end
    end
    
    out{ii}.peak_val = p_val;
    
    out{ii}.win = win;
    out{ii}.win_time = (-win:win)/mode(diff(behav.time));
end


%% for checking the process
% c_ord = linspecer(50);
% 
% figure(101)
% subplot(2,4, [1:4])
% cla
% hold on
% plot(wake_tvec, P_proj(ii,:));
% 
% plot(behav.time/1000, behav.position(:,1))
% 
% plot(out{ii}.loc_time/1000, out{ii}.peak_val, 'xk')
% 
% for ii = 1:length(out)
%     if ii >4; continue; end
%     subplot(2,4,4+ii)
%     cla; hold on
%     for ip = 1:length(out{ii}.loc)
%     plot((-win:win)/mode(diff(behav.time)), out{ii}.loc_mat(ip,:), 'color',[c_ord(ii,:) .5], 'linewidth', 2*(out{ii}.peak_val(ip)./max(out{ii}.peak_val)))
%     end
% 
% end

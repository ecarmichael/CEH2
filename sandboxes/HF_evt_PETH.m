%% sandbox_HF_OE_Peth





%% load the evts

evts_dir = ('/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/GoNoGo/HF2b2_D3/HF2b2_2025-12-18_13-47-24_D3_2_opto_only/Record Node 112/experiment1/recording1/events/Intan_RHD_USB-108.Rhythm Data/TTL') ;
csc_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/GoNoGo/HF2b2_D3/HF2b2_2025-12-18_13-47-24_D3_2_opto_only/Record Node 117'; 
phy_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/Kilo_inter/HF2b2_D3';
vr_fname = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/GoNoGo/HF2b2_D3/HF2b2_D3_2025-12-18_14-08-34.csv'; 

data = HF_preprocess(phy_dir, csc_dir, evts_dir, vr_fname, [1]); 

%% figure showing all the events

figure(1010)
clf
hold on
plot(data.S); 
plot(data.csc.tvec, (data.csc.data/100)-2, 'b')
plot(data.vr.pos.tvec, data.vr.pos.data/5, 'g', 'LineWidth',2)
plot(data.licks.tvec, (data.licks.data./5)-10, 'Color', [0.6 .6 .6])
ylim([-10 length(data.S.t)+1])
xline(data.evts.t{1}(:), 'm', 'Rewards')
xline(data.evts.t{2}(:), 'r')
xline(data.evts.t{3}(:), 'b')
xline(data.evts.t{4}(:), 'c', 'video')

xline(data.vr.evt.t{contains(data.vr.evt.label, 'Collision with Rwd1')}, 'c', 'LineStyle', '--')


figure(1012)
clf
hold on
plot(data.S); 
plot(data.csc.tvec, (data.csc.data/100)-2, 'b')
plot(data.vr.pos.tvec, data.vr.pos.data/5, 'g', 'LineWidth',2)
plot(data.licks.tvec, (data.licks.data./5)-10, 'Color', [0.6 .6 .6])
ylim([-10 length(data.S.t)+1])
xline(data.evts.t{1}(:), 'm', 'Rewards')
xline(data.evts.t{2}(:), 'r')
xline(data.evts.t{3}(:), 'b')
xline(data.evts.t{4}(:), 'c', 'video')

xline(data.vr.evt.t{contains(data.vr.evt.label, 'Collision with Rwd1')}, 'c', 'LineStyle', '--')

xlim([data.vr.pos.tvec(1) data.vr.evt.t{contains(data.vr.evt.label, 'Collision with Rwd1')}(1)+5])
%% peth per cell
c_red = [0.9153    0.2816    0.2878]; 


window = [-.250 .250]; 
bin_s = .025;
evt_t = data.evts.t{2} ;

% isolate events of a certain length
e_d = evt_t(2,:) - evt_t(1,:); 

ITIs = unique(round(e_d, 3)); 

for iTi = 1:length(ITIs)

    evt_t = data.evts.t{2} ;

k_idx = round(e_d,3) == ITIs(iTi) ;  

evt_t(:,~k_idx) = []; 
    % 
    % cfg = []; 
    % cfg.dt = 0.025; 
    % cfg.gauss_window = .01; 
    % cfg.gauss_sd = 0.0025;

    for iS = 1:length(data.S.t)

        % figure(iS+300)
        % clf

        S = data.S;
        S.t = [];
        S.t{1} = data.S.t{iS};
        S.label = [];
        S.label{1} =  data.S.label{iS};

        evt_t = sort(evt_t(1,:));

        cfg_peth = [];
        cfg_peth.window = [-.25 .25];
        cfg_peth. plot_type = 'raw';
        cfg_peth.dt = 0.001;
        cfg_peth.gauss_window = .025;
        cfg_peth.gauss_sd = .0025;
        cfg_peth.shuff = 500;
        cfg_peth.t_on = mode(e_d(k_idx));
        cfg_peth.rec_color = c_red;
        % cfg_peth.plot = 'off'; 
        [peth_S{iS,iTi}, peth_IT{iS,iTi},peth_gau{iS,iTi}, ~, ~, ~, ~, ~, ~,peth_T{iS, iTi}] = SpikePETH_Shuff(cfg_peth, S, evt_t);
    end
end

%% plot the PETHS together
reds = hot(10); 

for iS = 1:size(peth_gau,1)

    figure(iS);
    clf
    subplot(2,1,1)
    cla
    hold on;
    % Plot the PETH for each ITI
    offset = 0; 
    for iT = 1:size(peth_gau,2)
        u_val = unique(peth_T{iS, iT});
        for iV = 1:length(u_val)
            this_idx = peth_T{iS, iT} == u_val(iV);
            if isempty(this_idx)
                continue
            end
            plot(peth_S{iS, iT}(this_idx), peth_T{iS, iT}(this_idx)+0.5 + offset,'.', 'color', reds(iT,:), 'MarkerSize', 10)
            % disp(mode(peth_T{iS, iT}(this_idx)+ offset))
        end
        y_lim = ylim; 
        offset = y_lim(2); 

    end

    subplot(2,1,2)
        hold on;
    % Plot the mean activity for each ITI
    for iT = 1:size(peth_gau,2)
        plot(peth_IT{iS,iT}, mean(peth_gau{iS,iT},2, 'omitnan')', 'Color', reds(iT,:), 'LineWidth', 1.5);

    end
    xlabel('Time (s)');
    ylabel('PETH');
    title(['PETH for Cell ' num2str(iS) ' - ITI: ' num2str(ITIs(iT))]);
    legend({num2str(ITIs')}, 'Box', 'off')



end


% loop over events and remove the offset. 
% for ii = 1:length(evts.t)
%     evts.t{ii} = evts.t{ii} - csc.tvec(1); 
% end
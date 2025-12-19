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
window = [-1 5]; 
bin_s = .025;
evt_t = data.evts.t{2} ;

% isolate events of a certain length
e_d = evt_t(2,:) - evt_t(1,:); 

% k_idx = .08> e_d;   

% evt_t(:,~k_idx) = []; 

    cfg = []; 
    cfg.dt = 0.025; 
    % cfg.gauss_window = .01; 
    % cfg.gauss_sd = 0.0025;

for iS = 1:length(data.S.t)

    figure(iS+300)
    clf

    S = data.S;
    S.t = [];
    S.t{1} = data.S.t{iS};
    S.label = [];
    S.label{1} =  data.S.label{iS};

    evt_t = sort(evt_t(1,:)); 

 
        cfg_peth = [];
        cfg_peth.window = [-.5 .5];
        cfg_peth. plot_type = 'raw';
        cfg_peth.dt = 0.01;
        cfg_peth.gauss_sd = .025;
        cfg_peth.shuff = 500; 
        [~,~,this_peth] = SpikePETH_Shuff(cfg_peth, S, evt_t );
end


% loop over events and remove the offset. 
% for ii = 1:length(evts.t)
%     evts.t{ii} = evts.t{ii} - csc.tvec(1); 
% end
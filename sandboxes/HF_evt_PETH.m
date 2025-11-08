%% sandbox_HF_OE_Peth





%% load the evts

evts_dir = ('C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Wheel\test_data\SOM2_2025-10-22_15-43-05_opto_test_wide_probe4\Record Node 113\experiment1\recording1\events\Intan_RHD_USB-100.Rhythm Data\TTL') ;
csc_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Wheel\test_data\SOM2_2025-10-22_15-43-05_opto_test_wide_probe4\Record Node 118'; 
phy_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Wheel\Kilo_inter\SOM2_opto4';

data = HF_preprocess(phy_dir, csc_dir, evts_dir); 

%% peth per cell
window = [-1 5]; 
bin_s = .025;
evt_t = data.evts.t{2} ;

% isolate events of a certain length
e_d = evt_t(2,:) - evt_t(1,:); 

k_idx = .08> e_d;   

evt_t(:,~k_idx) = []; 

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
        cfg_peth.window = [-.25 .25];
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
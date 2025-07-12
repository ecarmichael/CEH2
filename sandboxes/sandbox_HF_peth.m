%% sandbox_HF_peth_summary





%% load some data

S = MS_Load_NTT([]);


evts = MS_LoadEvents();


%% split into recording blocks

start_idx = find(contains(evts.label, 'Starting Recording'));
stop_idx = find(contains(evts.label, 'Stopping Recording'));

Rec = []; 
for iR = length(evts.t{start_idx}):-1:1
    
   Rec{iR}.S = restrict(S, evts.t{start_idx}(iR),  evts.t{stop_idx}(iR)); 
      Rec{iR}.evts = restrict(evts, evts.t{start_idx}(iR),  evts.t{stop_idx}(iR)); 

    
    
end


%% loop over recordings and make a plot
c = {'r', 'b', 'm', 'k', 'y', 'g', 'c'};

for iR = 1:length(Rec)
    
%    HF_summary_peth(Rec{iR}.S, Rec{iR}.evts);
    cfg_S.spkColor = c{iR};
    PlotSpikeRaster2(cfg_S, Rec{iR}.S)
    
    
end
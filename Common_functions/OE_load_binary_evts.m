function [evts, timestamps, states, samples, words] = OE_load_binary_evts(evt_dir)
%% OE_load_binary_evts: loads the .npy binary TTL files [requires the npy-matlab toolbox]


p_dir = cd; 

cd(evt_dir)

timestamps = readNPY('timestamps.npy');
states = readNPY('states.npy');
samples = readNPY('sample_numbers.npy');
words = readNPY('full_words.npy');  % not used atm

%% use the states to determine the on and off of each state

state_id = unique(states); 
state_id(state_id <0) = []; 

evts = []; 
evts.type = 'ts'; 
for ii = 1:length(state_id)
    st_t = timestamps(states == state_id(ii)); 
    end_t = timestamps(states == -state_id(ii)); 


    evts.t{ii} = [st_t'; end_t']; 
    evts.label{ii} = num2str(state_id(ii)); 
end

cd(p_dir)
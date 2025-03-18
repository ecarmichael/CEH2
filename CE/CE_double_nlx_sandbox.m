function [csc1_out, csc_out2] = CE_double_nlx(Ca_dir1, Ca_dir2, nlx_dir, TTL)



if nargin < 4
    TTL = [5 3]; 
end


%% 

cd(nlx_dir)


cfg = []; 
cfg.fc = {'CSC1.ncs', 'CSC9.ncs', 'CSC12.ncs'};
cfg.desired_sampling_frequency = 2000; 

csc = MS_LoadCSC(cfg);

evts = MS_LoadEvents(); 

nlx1 = sort(unique([evts.t{TTL(1)} evts.t{TTL(1)+1}])); 
nlx2 = sort(unique([evts.t{TTL(2)} evts.t{TTL(2)+1}])); 


%% CA1
warning off
load([Ca_dir1 filesep 'ms.mat'])

ms_1 = MS_Ca_good_cells(ms); 
clear ms; 

% CA2
load([Ca_dir2 filesep 'ms.mat'])

ms_2 = MS_Ca_good_cells(ms); 
clear ms; 

warning on


%% restrict sleep section
ms_1_start_idx = [length(ms_1.tvecs{1}), length(ms_1.tvecs{1}) + length(ms_1.tvecs{2})]; 

ms1.time = ms_1.time(ms_1_start_idx(1)+1:ms_1_start_idx(2)); 
ms1.detrendRaw = ms_1.detrendRaw(ms_1_start_idx(1)+1:ms_1_start_idx(2),:); 
ms1.Binary = ms_1.Binary(ms_1_start_idx(1)+1:ms_1_start_idx(2),:); 


% trim to the sleep phase

ms_2_start_idx = [length(ms_2.tvecs{1}), length(ms_2.tvecs{1}) + length(ms_2.tvecs{2})]; 

ms2.time = ms_2.time(ms_2_start_idx(1)+1:ms_2_start_idx(2)); 
ms2.detrendRaw = ms_2.detrendRaw(ms_2_start_idx(1)+1:ms_2_start_idx(2),:); 
ms2.Binary = ms_2.Binary(ms_2_start_idx(1)+1:ms_2_start_idx(2),:); 

%% check the nunber of samples


fprintf('Ca1: %0d samples (%0.2fhrs)\n', length(ms1.time), (ms1.time(end) - ms1.time(1))/60/60);

fprintf('Ca2: %0d samples (%0.2fhrs)\n', length(ms2.time), (ms2.time(end) - ms2.time(1))/60/60);

fprintf('NLX %s: %0d samples (%0.2fhrs)\n','1', length(nlx1), (nlx1(end) - nlx1(1))/60/60);
fprintf('NLX %s: %0d samples (%0.2fhrs)\n','2', length(nlx2), (nlx2(end) - nlx2(1))/60/60);


for ii = 3:length(contains(evts.label, 'TTL'))
    
fprintf('NLX %s: %0d samples (%0.2fhrs)\n',evts.label{ii}, length(evts.t{ii}), (evts.t{ii}(end) - evts.t{ii}(1))/60/60);

end

%% 

if length(ms1.time)  < length(nlx1)
% remove the last TTL as is the case for   JKA_HPC_05\customEntValHere\2025_02_17  
    if diff(nlx1(end-1:end)) > mode(diff(nlx1))
        fprintf('last point has a jump (%.2fsec vs mode: %.2fsec). removing...\n', diff(nlx1(end-1:end)), mode(diff(nlx1)))
        nlx1(end) = []; 
    end
    
    
end
        

% more compicated version with dropped frames. 
% o_idx = diff(ms2.time) > mode(diff(ms2.time))*2; 

% if sum(o_idx) > 0
    
    
    
%    m2_interp = interp1(
    
    
%% output the MS1 intermediate data


MS_SWR_detector(csc, 'CSC1.ncs', 1)

%% proof of concept flash triggerd average



% this_ms = ms1; 


%% quick behaviour check

load([Ca_dir2 filesep 'behav_enc.mat'])


load([Ca_dir2 filesep 'behav_rec.mat'])

%%
ms2_enc.time = ms_2.time(1:length(ms_2.tvecs{1})); 
ms2_enc.Binary = ms_2.Binary(1:length(ms_2.tvecs{1}),:); 
ms2_enc.detrendRaw = ms_2.detrendRaw(1:length(ms_2.tvecs{1}),:); 
ms2_enc.deconv = ms_2.deconv(1:length(ms_2.tvecs{1}),:); 

behav_enc_a = MS_align_data(behav_enc, ms2_enc);

move_idx = behav_enc_a.speed > 5; 


for ii  =50:100
    figure(ii)
[rate_m, occ_mat ] = MS_decon_rate_map(ms2_enc.deconv(move_idx,ii),ms2_enc.time(move_idx),behav_enc_a.position(move_idx,:),2.5, 1, 1); 

end


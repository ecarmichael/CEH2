%%%%%%%%%%%%%%%%%%%%%%%%%% TUTORIAL ON CELL ASSEMBLY DETECTION %%%%%%%%%%%%%%%%%%%%%%%%%

% LOAD DATA
% load('test_data.mat');  

%% get the data into the spM format (spike times sorted with extra bins filled in with nans;

keep_idx = 1:(size(ms.deconv,2) - floor(size(ms.deconv,2)/3)); 
csp = ms.deconv(:,keep_idx)./ms.denoise(:,keep_idx); 
tvec = ms.time/1000; 

% get the max length of the number of 'spikes'
max_len = 0; 
for ii = size(csp,2):-1:1
   this_st{ii} = tvec((csp(:,ii) > 0.1)); 
   if length(this_st{ii})>max_len
       max_len = length(this_st{ii}); 
   end
   if isempty(this_st{ii})
       keep_idx(ii) = 0; 
   else
       keep_idx(ii) = 1; 
   end
end
% remove inactive cells
this_st(~keep_idx) = [];

% create a maxtrix that is nCells x nSpikeTimes with NaN padding. 
spM = NaN(length(this_st), max_len); 

for ii = size(this_st,2):-1:1
   spM(ii,1:length(this_st{ii})) = this_st{ii};  

end

%% run the CAD

nneu=size(spM,1);  % nneu is number of recorded units

BinSizes=[0.5 0.75 1 1.5  2 3 4 5];
MaxLags =[ 10  10  10 10  10 10  10 10];

% ASSEMBY DETECTION
[assembly]=Main_assemblies_detection(spM,MaxLags,BinSizes);

%% %%%%%%%%%%%%%%%%%%%%%%%% VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%

% ASSEMBLY REORDERING
[As_across_bins,As_across_bins_index]=assemblies_across_bins(assembly,BinSizes);

display='raw';
% display='clustered';

% VISUALIZATION
[Amatrix,Binvector,Unit_order,As_order]=assembly_assignment_matrix(As_across_bins, nneu, BinSizes, display);

%% %%%%%%%%%%%%%%%%%%%%%%%% PRUNING %%%%%%%%%%%%%%%%%%%%%%%%
clf
% PRUNING: criteria = 'biggest';
criteria = 'biggest';
[As_across_bins_pr,As_across_bins_index_pr]=pruning_across_bins(As_across_bins,As_across_bins_index,nneu,criteria);

display='raw';
[Amatrix,Binvector,Unit_order,As_order]=assembly_assignment_matrix(As_across_bins_pr, nneu, BinSizes, display);

%%
clf

% PRUNING: criteria = 'distance';
criteria = 'distance';
% th=0.7;
th=0.3;

style='pvalue';
% style='occ';

[As_across_bins_pr,As_across_bins_index_pr]=pruning_across_bins(As_across_bins,As_across_bins_index,nneu,criteria,th,style);
% display='raw';
display='clustered'
[Amatrix,Binvector,Unit_order,As_order]=assembly_assignment_matrix(As_across_bins_pr, nneu, BinSizes, display);


%% %%%%%%%%%%%%%%%%%%%%%%%% ASSEMBLY ACTVATION %%%%%%%%%%%%%%%%%%%%%%%%
clf

criteria = 'biggest';
[As_across_bins_pr,As_across_bins_index_pr]=pruning_across_bins(As_across_bins,As_across_bins_index,nneu,criteria);

lagChoice = 'beginning';
% lagChoice='duration';

act_count = 'full';
[assembly_activity]=Assembly_activity_function(As_across_bins_pr,assembly,spM,BinSizes,lagChoice,act_count);


m = ceil(length(assembly_activity)/2);
n = 2; 
for i=1:length(assembly_activity)
    subplot(m,n,i)
    plot(assembly_activity{i}(:,1),assembly_activity{i}(:,2));
    hold on
end














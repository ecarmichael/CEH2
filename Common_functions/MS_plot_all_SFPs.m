function hi = MS_plot_all_SFPs(SFPs_in)
%% MS_plot_all_SFPs:
%
%
%
%    Inputs: 
%    -  SPFs_in [X x Y x nSamples]  3d array with the spatial footprints of
%    the miniscope cells. should be the 'SPFs' subfield in the 'ms' structure from  
%
% % % %    - cutoff [double]  cut off value for normalized data.  Default is 0.
% % % %    put in a value between 0-1 to determine how much of the SFP remains.  
% % % %
%
%
%    Outputs: 
%    - hi : [handle]  handle for the imagesc. 
%
%
%
%
% EC 2020-10-26   initial version 
%
%
%
%% initialize

% if nargin <2
%     cutoff = 0; %default cutoff for converting to NaNs.  Can increase if you want to remove some of the lower values. 
% end
%% process the data
% % normalize the data, can use user input to determine cutoffs
% 
% SFPs_norm = SFPs_in./max(max(SFPs_in)); 
% 
% SFPs_norm(SFPs_norm <= cutoff) = NaN; 

% set background to match (needs to happen before the other nan "layers"
imagesc(zeros(size(SFPs_in,1), size(SFPs_in,2)))

for ii = 1:size(SFPs_in,3)
    hold on
    this_spf = SFPs_in(:,:,ii); 
    this_spf(this_spf==0) = NaN; 
        hi = imagesc(this_spf);
        set(hi,'alphadata',~isnan(this_spf));
end
axis xy
xlim([0 size(SFPs_in,2)])
ylim([0 size(SFPs_in,1)])


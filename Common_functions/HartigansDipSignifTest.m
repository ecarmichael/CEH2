function [dip, p_value, xlow,xup] = HartigansDipSignifTest(data_in, nboot, plot_flag)
%% wrapper function for Hartigan's dip test from Mechler / Nic Price using 
% the xpdf formating and significance tests from: 
% https://stackoverflow.com/questions/20815976/testing-for-unimodal-unimodality-or-bimodal-bimodality-distribution-in-matla and 
% https://snl.salk.edu/~jude/waveform_public/HartigansDipSignifTest.m

%% init


if nargin <2
    nboot = 500;
    plot_flag = 0; 
elseif nargin < 3
    plot_flag = 1; 
end


%% convert to xpdf

  x2 = reshape(data_in, 1, numel(data_in));
  [n, b] = hist(x2, 40);
  % This is definitely not probability density function
  xpdf = sort(x2);
  % downsampling to speed up computations
%   xpdf = interp1 (1:length(xpdf), xpdf, 1:1000:length(xpdf)); 
  
  %% run the test
  
[dip,xlow,xup, ~, ~, ~, ~, ~]=HartigansDipTest(xpdf);

  %% get the sig values from bootstrapping. 
  
  N=length(xpdf);

% calculate a bootstrap sample of size NBOOT of the dip statistic for a uniform pdf of sample size N (the same as empirical pdf)
boot_dip=[];
for i=1:nboot
   unifpdfboot=sort(unifrnd(0,1,1,N));
   [unif_dip]=HartigansDipTest(unifpdfboot);
   boot_dip=[boot_dip; unif_dip];
end
boot_dip=sort(boot_dip);
p_value=sum(dip<boot_dip)/nboot;

if plot_flag
    figure; clf;
    [hy,hx]=hist(boot_dip);
    bar(hx,hy,'k'); hold on;
    xline(dip, 'r:');
    
    
    
    
end
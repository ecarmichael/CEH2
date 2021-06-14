

%%  compute R2

win = 15;

for ii = 1:45
    this_bin = decoding.REM_decoded_probabilities(:,ii:ii+win);
%     for jj = 1:size(this_bin,2)
%         if isnan(sum(this_bin(:,jj))); continue; end
%         
%         this_bin(this_bin(:,jj) < max(this_bin(:,jj)),jj) = 0;
%     end
   P{ii}= polyfit(repmat(0:win, size(this_bin,1), 1), this_bin, 1);
end


%% 

col = randperm(size(this_bin,1)); 

for kk = size(this_bin,2):-1:1
    col = randperm(size(this_bin,1)); 

    this_shuff(:,kk) = this_bin(col,kk);
end

figure(101)
subplot(2,1,1)
imagesc(this_bin)

subplot(2,1,2)
imagesc(this_shuff)

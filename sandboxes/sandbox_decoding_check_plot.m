figure

imagesc(1:length(decoding.REM_decoded_position), min(decoding.REM_decoded_position):1:max(decoding.REM_decoded_position), decoding.REM_decoded_probabilities);
hold on
plot(1:length(decoding.REM_decoded_position), decoding.REM_decoded_position, '*', 'color', 'r')


for ii = 1:length(replay_idx)
    rectangle('position', [replay_idx(ii), min(decoding.REM_decoded_position)-1, 20, max(decoding.REM_decoded_position)-1], 'edgecolor', 'w')
end



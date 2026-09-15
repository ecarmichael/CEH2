function data = MS_h5_to_csv(fname, data_name); 


h_info = h5info(fname)

data = h5read(fname, '/df_with_missing/table');

data = data.values_block_0; 

writematrix(data, strrep(fname, 'h5', 'csv')); 

fprintf('Data from %s  saved back to %s\n', fname, strrep(fname, 'h5', 'csv'))
function data = MS_h5_to_csv(fname, data_name); 


h5info(fname)

data = h5read(fname, 'data_name');

writematrix(data, strrep(fname, 'h5', 'csv')); 

fprintf('Data from %s  saved back to %s\n', fname, strrep(fname, 'h5', 'csv'))
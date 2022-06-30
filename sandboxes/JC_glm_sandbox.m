%% JC GLM sandbox

load('GLM_matrix.mat')

%% convert to table

data_table  = array2table(GLM_firing_diff_matrix, 'VariableNames', {'isplace', 'MI', 'post_pre', 'post_task'});

data_table.isplace = logical(data_table.isplace);
col = jet(2);  
col_mat = zeros(length(data_table.isplace), 3); 
col_mat(data_table.isplace1:3) = col(1,:); 
summary(data_table) % print some basic stats

scatter3(data_table.post_task, data_table.post_pre, data_table.MI, [], data_table.isplace); 
xlabel('post_task');
ylabel('post_pre');
zlabel('MI')
legend({'Place'; 'non-place'});

%% make a simple GLM and compare it to a null

% build the generalized linear regression model. 

null_mdl = fitglm(data_table,('isplace ~ 1'),'link', 'logit', 'Distribution','binomial'); 

MI_mdl = fitglm(data_table,('isplace ~ 1 + MI'),'Distribution','binomial')
figure(101)
scatter(data_table.MI, data_table.post_pre, 10, data_table.isplace, 'filled')

% MI_only = fitglm(data_table,('isplace ~ MI'),'Distribution','binomial');
place_mdl = fitglm(data_table,('isplace ~ MI * post_task * post_pre'),'Distribution','binomial')
% anova(place_all)
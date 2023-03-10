%% JC GLM sandbox

load('GLM_matrix.mat')

%% convert to table

data_table  = array2table(GLM_firing_diff_matrix, 'VariableNames', {'isplace', 'MI', 'post_pre', 'post_task'});

data_table.isplace = logical(data_table.isplace);
col =[0.9153 0.2816 0.2878; 0.3467 0.5360 0.6907] ;  
col_mat = zeros(length(data_table.isplace), 3); 
for ii  = length(data_table.isplace):-1:1
    if data_table.isplace(ii)
        col_mat(ii,1:3) = col(1,:);
    else
        col_mat(ii,1:3) = col(2,:);
    end
end
summary(data_table) % print some basic stats


%% make a simple GLM and compare it to a null

% build the generalized linear regression model. 
% null model 
null_mdl = fitglme(data_table,('isplace ~ 1'),'link', 'logit', 'Distribution','binomial'); 

% MI only 
MI_mdl = fitglme(data_table,('isplace ~ 1 + MI'),'Distribution','binomial')

% all factors
place_mdl = fitglme(data_table,('isplace ~ MI * post_task * post_pre'),'Distribution','binomial')

% visualize 
figure(101)
clf
subplot(1,2,1)
hold on
scatter(data_table.post_pre(data_table.isplace),data_table.MI(data_table.isplace), 50,col_mat(data_table.isplace,:), 'filled')
scatter(data_table.post_pre(~data_table.isplace),data_table.MI(~data_table.isplace), 50,col_mat(~data_table.isplace,:))
ylabel('MI')
xlabel('post_pre');
legend({'Place'; 'non-place'});

subplot(1,2,2)
hold on
scatter3(data_table.post_task(data_table.isplace), data_table.post_pre(data_table.isplace), data_table.MI(data_table.isplace), 50, col_mat(data_table.isplace,:), 'filled'); 
scatter3(data_table.post_task(~data_table.isplace), data_table.post_pre(~data_table.isplace), data_table.MI(~data_table.isplace), 50, col_mat(~data_table.isplace,:)); 
view(-45,25)
xlabel('post_task');
ylabel('post_pre');
zlabel('MI')
legend({'Place'; 'non-place'}, 'Location', 'northeast');




%% compare the place to the null models using AIC/BIC (smaller number is better)
% compare(MI_mdl, null_mdl)
% compare(place_mdl, null_mdl)
compare(place_mdl, MI_mdl)

% since MI only model has lower AIC/BIC is it the 'better' model. 
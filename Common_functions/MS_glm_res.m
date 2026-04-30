function [res, mdl, tbl] = MS_glm_res(data1, data2)
%% MS_glm_res:  simple function to gresses out data2 from data1 and output the residual.  assumes more is data1 ~ data2

% convert to table
tbl = table(data1, data2, 'VariableNames',{'data1', 'data2'}); 

% Fit a linear model to data1 using data2 as predictor
mdl = fitglm(tbl, 'data1 ~ data2'); 

% get the residual
res = mdl.Residuals.Raw;

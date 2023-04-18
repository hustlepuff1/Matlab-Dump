%% Linear Regression (Simplest Implementation)
%% Generating Data
%%
% Taking any random data
X=randn(100,1);

% Adding noise to equation y = 3x + 7
Y=3*(X + 0.3*randn(100,1)) + (7 + 0.3*randn(100,1));

% If data is not normalized, Zscore Normalization (Advised for SGD)
% data=[zscore(X),Y];
data=[X,Y];

% Define a symbolic varible for function plot
syms x
%% CV Partition
% * 70% to training set + 30% to testing set 

[train_set,test_set ] = holdout(data,70 );

X_train=train_set(:,1:end-1); Y_train=train_set(:,end);
X_test=test_set(:,1:end-1); Y_test=test_set(:,end);

% Number of training instances
N=length(X_train)

% Number of testing instances
M=length(X_test)
%% 1.  Using Direct Method
% * Append a vectors of one to _X_train _for calculating *bias*.
%%
W=pinv([ones(N,1) X_train])*Y_train
%% Mean Square Error

predicted_values=[ones(M,1) X_test]*W;
mse1=sqrt(mean((predicted_values-Y_test).^2))
%% Plot

figure
hold on
scatter(X_test,Y_test)
fplot(W(1)+W(2)*x)
xlabel({'X_1'})
ylabel({'Y'})
title({'Regression Using Direct Method'})
xlim([-3 3])
hold off
%% 2. Using Inbuilt MATLAB Function
%%
Test_mdl = fitlm(X_train,Y_train);
W=Test_mdl.Coefficients{:,1}
%% Mean Square Error

predicted_values=predict(Test_mdl,X_test);
mse2=sqrt(mean((predicted_values-Y_test).^2))
%% Plot

figure
hold on
scatter(X_test,Y_test)
fplot(W(1)+W(2)*x)
xlabel({'X_1'});
ylabel({'Y'});
title({'Regression Using Inbuilt MATLAB Function'});
xlim([-3 3])
hold off
%% 3. Using Stochastic Gradient Descent
% * Learning Parameter, alpha = 0.1
% * Append a vectors of one to _X_train _for calculating *bias*.
% * Tolerence = 10^-5
%%
X_train=[ones(N,1), X_train];

W=zeros(size(X_train,2),1);
W_old=ones(size(X_train,2),1);

while(norm(W_old-W) > 10^-5)
    W_old=W;
    W = W - 0.1/N*X_train'*(X_train*W - Y_train);
end
W
%% Mean Square Error

predicted_values=[ones(length(X_test),1),X_test]*W;
mse3=sqrt(mean((predicted_values-Y_test).^2))
%% Plot

figure
hold on
scatter(X_test,Y_test)
fplot(W(1)+W(2)*x)
xlabel({'X_1'});
ylabel({'Y'});
title({'Regression using Stochastic Gradient Descent'});
xlim([-3 3])
hold off
%% Comparing Mean Square
%%
figure
mse=[mse1,mse2,mse3];
plot(mse,'--^m')
title('Comparison: Mean Square Error')
ylabel({'MSE \rightarrow'});
xlabel({'Different methods   [--Direct_{Method} --Inbuilt_{Function} --SGD_{Method}]'});
xlim([0, 4]); ylim([min(mse)-1e-6, max(mse)+1e-6]);
%%
% Created by Bhartendu, Machine Learning & Computing

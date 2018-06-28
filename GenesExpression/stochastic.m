function [ X_stoch ] = stochastic( X )

D = diag(sum(X,2));
X_stoch = D\X;

end


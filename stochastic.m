function [ X_stoch ] = stochastic( X )
% Converts a given matrix to a right stochastic matrix.
% 
% Input:
%     X - real matrix with nonnegative entries  
% 
% Output:
%     X_stoch - right stochastic martix, with each row summing to 1  
%--------------------------------------------------------------------------

D = diag(sum(X,2));
X_stoch = D\X;

end


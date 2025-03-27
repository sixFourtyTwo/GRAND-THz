% H: m x n; m >= n
% b: m x 1
% When m > n, there is a residual in the solution of ||Hx-b||_2
function [Q, R] = qrd1(H)

[m,n] = size(H);
R = complex(zeros(min(m,n),n),zeros(min(m,n),n));
Q = H;

for i=1:n % 1 <= pv <= n  
    % orthognalize Q,R at col pos i
    R(i,i) = sqrt(Q(:,i)'*Q(:,i));
    Q(:,i) = Q(:,i)/R(i,i);    
    for j=i+1:n
        R(i,j) = Q(:,i)'*Q(:,j);        
        Q(:,j) = Q(:,j) - R(i,j)*Q(:,i);        
    end
end

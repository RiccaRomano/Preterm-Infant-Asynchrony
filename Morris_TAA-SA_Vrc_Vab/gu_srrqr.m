function [cols] = gu_srrqr(A, k)
%% function [cols] = qu_srrqr(A, k)
%
% This function determines a good set of k "maximally independent"
% columns of A by the Gu-Eisenstat strong rank revealing
% QR factorization.
%
% Input:
%   A, an m by n real matrix
%   k, the number of columns to choose
%
% Output:
%   cols, the indices of k columns of A
%
%  The submatrix of good columns is A(:, cols) and it satifies the
%  following:
%
%  Let A1 = A(:,cols), A2 = A(:, setdiff(1:n, cols), 
%  C = A1*pinv(A1)*A2-A2, [s] = svd(A); s1 = svd(A1); s2 = svd(C), 
%  then
%  for i = 1,...,k, s1(i) >= s(i)/sqrt(1+k(n-k))
%  also
%  abs(pinv(A1)*A2) has no entry > 1
%  and
%  for j = 1,..., n-k
%  s2(j) <= s(k+j)*sqrt(1+k(n-k))
%    
%%  


        
% Initializations
n = size(A, 2);
r = A; 
%P will keep track of column permutations
[q, r, P] = qr(r, 0);
increasefound = true;
while (increasefound)
        [q, r] = qr(r, 0); %triangularize r
        %partition r into
        %  k (n-k)
        %| a   b |k
        %| 0   C |(m-k)
        a = r(1:k,1:k); 
        ainvb = a\r(1:k,k+1:end); %form a^{-1}b
        C = r(k+1:end, k+1:end);
        
        %compute the column norms of C
        gamma = zeros(size(C, 2), 1);
        for ccol = 1:size(C, 2)
            gamma(ccol) = norm(C(:,ccol), 2);
        end
        
        %find row norms of a^-1
        [l, u] = lu(a');
        w = zeros(k, 1);
        for arow = 1:k
            tmp = zeros(k, 1);
            tmp(arow) = 1;
            w(arow) = norm(u\(l\tmp), 2);
        end
        
        %find indices i and j that maximize 
        %ainv(i,j)^2 + (w(i)*gamma(j))^2
        tmp = w*gamma';
        fij = ainvb.^2 + tmp.^2;
        [f, i, j] = maxmat(abs(fij));
        
        %sqrt(f) gives the factor increase in |det(a)| 
        %when columns i and j are swapped
        if (f > 1)
            %we can increase |det(a)|
            increasefound = true;
            r(:,[i j+k]) = r(:, [j+k i]);  % permute columns i and j
            P([i j+k]) = P([j+k i]);
        else
            %we're done
            increasefound = false;
        end
end
    
cols = P(1:k);

    
function [m, i, j] = maxmat(A)
% find a maximum element of the matrix A and return
% it's row and column index
  k = size(A, 1);
  [m, ind] = max(A(:)); %A(:) puts the elements of A in a column vector
  i = mod(ind-1, k) + 1;
  j = floor((ind-1)/k) + 1;
function [p,Perm,err,Am] = cholesky_pivoted(A,etol)

%etol is the accepted error
%A is the input 2D matrix

%first I sort the A matrix so that diagonal matrix elements are in descending order 
dA = diag(A);
n = length(dA);
[SdA,order] = sort(dA,'descend');
Perm = eye(n);
Perm = Perm(:,order);

Ap = Perm'*A*Perm;

Am = zeros(size(A));
err = 1;
m = 1;
while err > etol
    
    l(m,m) = sqrt(SdA(m));
    if m == 1
        l(m,m+1:n) = Ap(m,m+1:n)/l(m,m);
        SdA(m+1:n) = SdA(m+1:n)-l(m,m+1:n).*l(m,m+1:n);
    end
    
end

for i = 1:m
    Am = Am + l(i,:)*l(i,:)';
end

end

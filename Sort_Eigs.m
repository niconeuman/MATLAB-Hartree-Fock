function [eigenvecs, eigenvals] = Sort_Eigs(C,epsilon)

nb = size(C,1);
[sorted_eigvals, index] = sort(diag(epsilon));
sorted_eigvecs = zeros(size(C));
    for n = 1:nb
        sorted_eigvecs(:,n) = C(:,index(n));
    end
    eigenvals = sorted_eigvals;
    eigenvecs = sorted_eigvecs;

end

function [p,PP,err,err2] = cholesky_four_ind(Inp_mat,etol)

%Input_matrix: is a four index electron repulsion matrix. Eventually it
%should not be necessary to calculate it. The cholesky vectors should be
%calculated from the four index matrix elements calculated on the fly.

%This is an approximated decomposition of a matrix (I start with two
%indices, then I will try to expand to 4)
%nvec: number of Cholesky vectors
%etol: error in the approximation
if size(Inp_mat,1) == size(Inp_mat,2)
    Dim = size(Inp_mat,1);
end

K = 1;

%chol_vec = zeros(Dim,Dim);
err = 10;
dA = diag(Inp_mat);
pi_ = (1:Dim);
l_i = cell(Dim,1);
    while (err > etol)
    
        [~,pivot] = max(dA(pi_(K:Dim))); %Find the maximum value and index of the diagonal of the remainder of the input_matrix
        tmp = pi_(pivot);
        pi_(pivot) = pi_(K);
        pi_(K) = tmp;
        l_i{K} = zeros(1,Dim);
        l_i{K}(pi_(K)) = sqrt(dA(pi_(K)));
        
       
        if K == 1
           for i = K + 1:Dim
                l_i{K}(pi_(i)) = 1/l_i{K}(pi_(K))*(Inp_mat(pi_(K),pi_(i)));
                dA(pi_(i))=dA(pi_(i))-l_i{K}(pi_(i))*l_i{K}(pi_(i));
           end
        
        else
            for i = K + 1:Dim
                sum_l_iK = 0;
                    for m = 1:K-1
                        sum_l_iK = sum_l_iK + l_i{m}(pi_(K))*l_i{m}(pi_(i));
                    end
                l_i{K}(pi_(i)) = 1/l_i{K}(pi_(K))*(Inp_mat(pi_(K),pi_(i))-sum_l_iK);
                dA(pi_(i))=dA(pi_(i))-l_i{K}(pi_(i))*l_i{K}(pi_(i));
            end
        end
        dA
        
        err = 0;
        for i = K + 1:Dim
            err = err + abs(dA(pi_(i)));
        end
        
        err2 = sum(dA(pi_(K+1:Dim)));
        K = K + 1;

    end
    p = cell2mat(l_i);
    
    PP = eye(Dim); PP = PP(:,pi_);
    %P = cell2mat(p{1:K});
end


function [nvec,chol_vec] = cholesky_two_ind(Inp_mat,nvec,~)

%Input_matrix: is a four index electron repulsion matrix. Eventually it
%should not be necessary to calculate it. The cholesky vectors should be
%calculated from the four index matrix elements calculated on the fly.

%This is an exact decomposition of a 2D matrix

if size(Inp_mat,1) == size(Inp_mat,2)
    Dim = size(Inp_mat,1);
end


chol_vec = zeros(Dim,Dim);

for k = 1:Dim

    if k == 1
        chol_vec(k,k) = Inp_mat(k,k)^0.5;
        for l = k+1:Dim
        chol_vec(l,k) = Inp_mat(l,k)/chol_vec(k,k);
        end
        
    else
    
        sum_hsq = 0;
        for j = 1:k-1
            sum_hsq = sum_hsq + chol_vec(k,j)^2 ;
        end
    
            chol_vec(k,k) = (Inp_mat(k,k) - sum_hsq)^0.5;
    
        for l = k+1:Dim
            sum_h_h = 0;
            for m = 1:k-1
                sum_h_h = sum_h_h + chol_vec(l,m)*chol_vec(k,m);
            end
            chol_vec(l,k) = 1/chol_vec(k,k)*(Inp_mat(l,k)-sum_h_h);
        end
    end
    
end    
end


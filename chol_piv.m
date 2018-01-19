function [L,err,I,J,M] = chol_piv(A,etol)

%A is a square positive semidefinite matrix
err = 1;
J = 1;
dA = diag(A);
M = length(dA);
%This programs sort of works. It outputs a matrix L which gives
%L*L' ~ A
%but I think it currently only works when the diagonal elements already are
%in descending order
while err > etol
    dA = diag(A);
    if J == 1
    [~,Jmax] = max(dA); %find maximum diagonal element
        if J == Jmax
            L(J,J) = sqrt(dA(J));
        else
            dA([J Jmax]) = dA([Jmax J]); %this permutes the elements
            L(J,J) = sqrt(dA(J));
        end
        for I = J+1:M
            L(I,J) = 1/L(J,J)*(A(I,J));
        end
    else
        L(J,J) = sqrt(A(J,J)-sum(L(J,1:J-1).^2));
        for I = J+1:M
            L(I,J) = 1/L(J,J)*(A(I,J)-sum(L(I,1:J-1).*L(J,1:J-1)));
        end
    end
    
    if J == M
        break
    end
    
    err = sum(abs(dA));
    %err = L(J,J)^2;
    
    
    J  = J+1;
    
end

end
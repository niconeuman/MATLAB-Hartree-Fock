function [MinE,E,ncycle,D,G,epsilon,F,Fprime] = SCF(H0,gabcd,S,Nel)
maxcycles = 30;
converged = 1;
ncycle = 0;

[V,D] = eig(S);
X = V*D^(-0.5)*V'; %transforms to an orthonormal basis

F = H0;
Fprime = X'*F*X;

[Cprime,epsilon] = eig(Fprime);

C = X*Cprime;
[C,epsilon] = Sort_Eigs(C,epsilon);

D = Build_Density(C,Nel);

E = zeros(maxcycles,1);

while ncycle < maxcycles && converged ~= 0
    ncycle = ncycle + 1;
    G = Build_Coulomb_Exchange(D,gabcd);
    F = H0 +  G;
    Fprime = X'*F*X;
    [Cprime,epsilon] = eig(Fprime);
    C = X*Cprime;
    [C,epsilon] = Sort_Eigs(C,epsilon);
    
    D = Build_Density(C,Nel);
    
    E(ncycle) = Fock_Energy(D,H0,F);
        if(ncycle > 1) && abs(E(ncycle)-E(ncycle-1)) < 1e-6
            converged = 0;
        end
end    
E = E(1:ncycle);
MinE = E(ncycle);
end
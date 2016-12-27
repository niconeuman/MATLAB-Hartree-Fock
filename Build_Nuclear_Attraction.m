function [VAB,nAtom,nbasis] = Build_Nuclear_Attraction(basis,AL,Z)
nAtom = size(Z,2);
nbasis = size(basis,2);
VAB = zeros(nbasis,nbasis);

for N = 1:nAtom
    for a = 1:nbasis
        for b = 1:nbasis
            for na  = 1:basis{a}.n
                for nb = 1:basis{b}.n
                    VAB(a,b) = VAB(a,b) + coulombg(basis{a}.g(na),basis{b}.g(nb),AL(N,1),AL(N,2),AL(N,3),Z(N))*basis{a}.c(na)*basis{b}.c(nb);
                end
            end
        end
    end
end

end
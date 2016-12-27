function D = Build_Density(eigenvecs,N)

nb = size(eigenvecs,1);
D = zeros(nb,nb);

for n = 1:nb
    for m = 1:nb
        for i = 1:N/2
            D(n,m) = D(n,m)+eigenvecs(n,i)*eigenvecs(m,i);
        end
    end
end

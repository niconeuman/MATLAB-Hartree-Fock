function G = Build_Coulomb_Exchange(D,gabcd)
nb = size(D,1);
G = zeros(nb,nb);

for a = 1:nb
    for b = 1:nb
        for c = 1:nb
            for d = 1:nb
                G(a,b) = G(a,b) + 0.5*D(c,d)*(2*gabcd(a,b,c,d)-gabcd(a,d,c,b));
            end
        end
    end
end

end
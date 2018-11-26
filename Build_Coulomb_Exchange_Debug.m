function [Gslow,JpK,JpK1,JpK2,Gfast] = Build_Coulomb_Exchange_Debug(D,gabcd)
nb = size(D,1);
Gslow = zeros(nb,nb);
G = zeros(nb,nb);

%Gslow
tic;

for i = 1:nb
    for j = 1:nb
        for k = 1:nb
            for l = 1:nb
                Gslow(i,j) = Gslow(i,j) + 0.5*D(k,l)*(2*gabcd(i,j,k,l)-gabcd(i,l,k,j));
            end
        end
    end
end
toc;

J = zeros(nb,nb);
K = zeros(nb,nb);
J1 = zeros(nb,nb);
K1 = zeros(nb,nb);
J2 = zeros(nb,nb);
K2 = zeros(nb,nb);

%JpK
tic;
for i = 1:nb
    for j = 1:nb
        for k = 1:nb
            for l = 1:nb
                J(i,j) = J(i,j) + D(k,l)*gabcd(i,j,k,l);
                K(i,l) = K(i,l) - 0.5*D(k,j)*gabcd(i,j,k,l);
            end
        end
    end
end
toc;
%JpK1
tic;
for i = 1:nb
    for j = 1:i
        for k = 1:nb
            for l = 1:nb
                J1(i,j) = J1(i,j) + D(k,l)*gabcd(i,j,k,l);
            end
        end
    end
end

for i = 1:nb
    for j = 1:nb
        for k = 1:nb
            for l = 1:i
                K1(i,l) = K1(i,l) - 0.5*D(k,j)*gabcd(i,j,k,l);
            end
        end
%        J(j,i) = J(i,j);
%        K(j,i) = K(i,j);
    end
end
toc;

%Jpk2
tic;
for i = 1:nb
    for j = 1:i
        for k = 1:nb
            for l = 1:k-1
                J2(i,j) = J2(i,j) + 2*D(k,l)*gabcd(i,j,k,l);
            end
            l = k;
                J2(i,j) = J2(i,j) + D(k,l)*gabcd(i,j,k,l);
        end
    end
end

%Just permutting the j and l loops, with their respective limits, does not change the result
for i = 1:nb
    for l = 1:i
        for k = 1:nb
            for j = 1:nb
                K2(i,l) = K2(i,l) - 0.5*D(k,j)*gabcd(i,j,k,l);
            end
        end
    end
end
toc;

%All possible permutations of indices (I switch the coulomb integral, and other indexes accordingly)
%Gslow(ij) = Gslow(ij) + 0.5*D(kl)*(2*gabcd(ijkl)-gabcd(ilkj));
%Gslow(ji) = Gslow(ji) + 0.5*D(kl)*(2*gabcd(jikl)-gabcd(jlki));
%Gslow(ij) = Gslow(ij) + 0.5*D(lk)*(2*gabcd(ijlk)-gabcd(iklj));
%Gslow(ji) = Gslow(ji) + 0.5*D(lk)*(2*gabcd(jilk)-gabcd(jkli));
%Gslow(kl) = Gslow(kl) + 0.5*D(ij)*(2*gabcd(klij)-gabcd(kjil));
%Gslow(kl) = Gslow(kl) + 0.5*D(ji)*(2*gabcd(klji)-gabcd(kijl));
%Gslow(lk) = Gslow(lk) + 0.5*D(ij)*(2*gabcd(lkij)-gabcd(ljik));
%Gslow(lk) = Gslow(lk) + 0.5*D(ji)*(2*gabcd(lkji)-gabcd(lijk));

%Then I order them in the indexes for G

%Gslow(ij) = Gslow(ij) + 0.5*D(kl)*(2*gabcd(ijkl)-gabcd(ilkj));
%Gslow(ij) = Gslow(ij) + 0.5*D(lk)*(2*gabcd(ijlk)-gabcd(iklj));
%Gslow(ji) = Gslow(ji) + 0.5*D(lk)*(2*gabcd(jilk)-gabcd(jkli));
%Gslow(ji) = Gslow(ji) + 0.5*D(kl)*(2*gabcd(jikl)-gabcd(jlki));
%Gslow(kl) = Gslow(kl) + 0.5*D(ij)*(2*gabcd(klij)-gabcd(kjil));
%Gslow(kl) = Gslow(kl) + 0.5*D(ji)*(2*gabcd(klji)-gabcd(kijl));
%Gslow(lk) = Gslow(lk) + 0.5*D(ij)*(2*gabcd(lkij)-gabcd(ljik));
%Gslow(lk) = Gslow(lk) + 0.5*D(ji)*(2*gabcd(lkji)-gabcd(lijk));

Jfactor = 2;
Kfactor = 0.5;
%Gfast
tic;
for i = 1:nb
    for j = 1:i
        if i == j
            factorij = 1;
        else
            factorij = 2;
        end

        for k = 1:i
            if (k == i)
                for l = 1:j
                    if l == j
                        factorlj = 1;
                    else
                        factorlj = 2;
                    end
                    if l == k
                        factorlk = 1;
                    else
                        factorlk = 2;
                    end
                    gijkl = gabcd(i,j,k,l);

                    G(i,j) = G(i,j) +  Jfactor*factorij*D(k,l)*gijkl;
                    G(k,l) = G(k,l) +  Jfactor*factorlk*D(i,j)*gijkl;
                    G(i,k) = G(i,k) - Kfactor*D(j,l)*gijkl;
                    G(j,l) = G(j,l) - Kfactor*factorlj*D(i,k)*gijkl;
                    G(i,l) = G(i,l) - Kfactor*D(j,k)*gijkl;
                    G(j,k) = G(j,k) - Kfactor*D(i,l)*gijkl;
                end
            else
                for l = 1:k
                    if l == j
                        factorlj = 1;
                    else
                        factorlj = 2;
                    end
                    if l == k
                        factorlk = 1;
                    else
                        factorlk = 2;
                    end
                    gijkl = gabcd(i,j,k,l);
                    G(i,j) = G(i,j) +   Jfactor*factorij*D(k,l)*gijkl;
                    G(k,l) = G(k,l) +   Jfactor*factorlk*D(i,j)*gijkl;
                    G(i,k) = G(i,k) - 2*Kfactor*D(j,l)*gijkl;
                    G(j,l) = G(j,l) - Kfactor*factorlj*D(i,k)*gijkl;
                    G(i,l) = G(i,l) - Kfactor*D(j,k)*gijkl;
                    G(j,k) = G(j,k) - Kfactor*D(i,l)*gijkl;
                end
            end
        end

    end
end
toc;
%for i = 1:nb
%    for j = 1:i
%        G(j,i) = G(i,j);
%    end
%end

Gfast = G;

JpK = J + K;
JpK1 = J1 + K1;
JpK2 = J2 + K2;
end

function GProd_Data = Plot_GProd(basis,Pcell,pcell,S,aa,bb)

%Aug 27th 2017
%This function takes to basis functions basis_a and basis_b, Pcell and
%pcell and plots all gaussian products between the primitives. It also
%plots the approximated integrals from Pcell, and pcell
%Maybe it needs more parameters, such as S

basis_a = basis{aa};
basis_b = basis{bb};
zz = linspace(-15,15,500)';
Gzz_sep = zeros(length(zz),basis_a.n*basis_b.n);
Gzz_accum = zeros(length(zz),1);
t = 1;
GProd_Data = zeros(basis_a.n*basis_b.n,7);
figure


    for nba = 1:basis_a.n
        g1 = basis_a.g(nba);
            for nbb = 1:basis_b.n
                g2 = basis_b.g(nbb);
                a = g1.alpha;
                b = g2.alpha;
                %Gaussian coefficient of the product
                p = a + b;
                q = a*b/p;
                RAB = [g1.x0-g2.x0;g1.y0-g2.y0;g1.z0-g2.z0];
                
                KAB = exp(-q*(RAB(1)^2+RAB(2)^2+RAB(3)^2));
                %Coordinates of the product gaussian
                Px = (g1.alpha*g1.x0 + g2.alpha*g2.x0)/p;
                Py = (g1.alpha*g1.y0 + g2.alpha*g2.y0)/p;
                Pz = (g1.alpha*g1.z0 + g2.alpha*g2.z0)/p;
                
                GProd_Data(t,:) = [a,b,p,Px,Py,Pz,KAB];
                
                [x, y, z] = ellipsoid(Px,Py,Pz,1/p,1/p,1/p,30);
                s = surf(x, y, z); hold on
                s.EdgeColor = 'none';
                s.FaceColor = 'r';
                s.FaceAlpha = 0.1*KAB;
                axis equal
                
                
                Gzz = KAB*exp(-p*(zz-Pz).^2);
                Gzz_sep(:,t) = Gzz;
                Gzz_accum = Gzz_accum + Gzz*basis_a.c(nba)*basis_b.c(nbb)*basis_a.g(nba).N*basis_b.g(nbb).N;
                
                t = t+1;
            end        
    end
    Gaprox = 1/(sqrt(pi/pcell{aa,bb})^3)*S(aa,bb)*exp(-pcell{aa,bb}*(zz-Pcell{aa,bb}(3)).^2);
    
    figure
    plot(zz,Gzz_sep,'-k'); hold on
    plot(zz,Gzz_accum,'-k','LineWidth',2); hold on
    plot(zz,Gaprox,'-r','LineWidth',2);
end
function T = gkinetic(g1,g2)

a = g1.alpha;
b = g2.alpha;
p = a+b;
%P = [a*g1.x0 + b*g2.x0,a*g1.y0 + b*g2.y0,a*g1.z0 + b*g2.z0]/p;
R12sq = (g1.x0-g2.x0)^2+(g1.y0-g2.y0)^2+(g1.z0-g2.z0)^2;

%RPA2 = sum((A-P).^2);

    Ex = gprod_1D(g1.x0,g1.alpha,g2.x0,g2.alpha);
    Ey = gprod_1D(g1.y0,g1.alpha,g2.y0,g2.alpha);
    Ez = gprod_1D(g1.z0,g1.alpha,g2.z0,g2.alpha);
    
    S = Ex*Ey*Ez*sqrt(pi/(g1.alpha + g2.alpha))^3*g1.N*g2.N;
    T = S*(a*b/p*(3-2*a*b/p*R12sq)); 


end
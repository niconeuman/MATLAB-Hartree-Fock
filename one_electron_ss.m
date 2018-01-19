function [S,T,P,p] = one_electron_ss(g1,g2,RAB,Px,Py,Pz)
%November 9th 2017

a = g1.alpha;
b = g2.alpha;
p = a+b;
q = a*b/p;
%P = [a*g1.x0 + b*g2.x0,a*g1.y0 + b*g2.y0,a*g1.z0 + b*g2.z0]/p;
R12sq = (g1.x0-g2.x0)^2+(g1.y0-g2.y0)^2+(g1.z0-g2.z0)^2;

   
    KAB = exp(-q*(RAB(1)^2+RAB(2)^2+RAB(3)^2));
    
    S = KAB*sqrt(pi/(g1.alpha + g2.alpha))^3;
    T = S*(a*b/p*(3-2*a*b/p*R12sq));
    P = [Px;Py;Pz];
end
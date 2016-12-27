function V = coulombg(g1,g2,Ax,Ay,Az,Z)

a = g1.alpha;
b = g2.alpha;
p = a+b;
P = [a*g1.x0 + b*g2.x0,a*g1.y0 + b*g2.y0,a*g1.z0 + b*g2.z0]/p;
A = [Ax,Ay,Az];

RPA2 = sum((A-P).^2);

Ex = gprod_1D(g1.x0,g1.alpha,g2.x0,g2.alpha);
Ey = gprod_1D(g1.y0,g1.alpha,g2.y0,g2.alpha);
Ez = gprod_1D(g1.z0,g1.alpha,g2.z0,g2.alpha);

V = -1*Z*2*pi/p*Boys(0,p*RPA2)*g1.N*g2.N;



end
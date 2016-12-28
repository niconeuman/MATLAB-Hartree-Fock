function orbital = Build_SOrbital(Ax,Ay,Az,alpha)
%This calculates the normalization constant for a gaussian orbital.
%I need to extend this for arbitrary angular momentum


tmp.x0 = Ax;
tmp.y0 = Ay;
tmp.z0 = Az;
tmp.alpha = alpha;
tmp.N = (2*alpha/pi)^.75;
tmp.L = 0;
orbital = tmp;


end
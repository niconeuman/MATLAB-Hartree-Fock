function orbital = Build_POrbital(Ax,Ay,Az,alpha)
%mx = 1, or my = 1 or mz = 1
tmp.x0 = Ax;
tmp.y0 = Ay;
tmp.z0 = Az;
tmp.alpha = alpha;
%L = 1;
tmp.N = 2^(1)*(2/pi)^.75*alpha^(0.75+0.5*(1))/(double_factorial(1,0,0))^.5;
tmp.L = 1;
orbital = tmp;


end
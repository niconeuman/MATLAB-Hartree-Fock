function orbital = Build_POrbital(Ax,Ay,Az,alpha)
%mx = 1, or my = 1 or mz = 1
tmp.x0 = Ax;
tmp.y0 = Ay;
tmp.z0 = Az;
tmp.alpha = alpha;
tmp.N = (2/pi)^.75*(2^1*alpha^((2+3)/4)/(double_factorial(1,0,0))^.5);
tmp.L = 1;
orbital = tmp;


end
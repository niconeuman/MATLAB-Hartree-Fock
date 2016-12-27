function orbital = Build_SOrbital(Ax,Ay,Az,alpha)

tmp.x0 = Ax;
tmp.y0 = Ay;
tmp.z0 = Az;
tmp.alpha = alpha;
tmp.N = (2*alpha/pi)^.75;
orbital = tmp;


end
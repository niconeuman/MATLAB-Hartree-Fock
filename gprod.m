function [KAB,p,q] = gprod(x1,y1,z1,alpha1,x2,y2,z2,alpha2)
p = alpha1+alpha2;
q = (alpha1*alpha2)/p;
%R1 = [x1;y1;z1];
%R2 = [x2;y2;z2];
%P = (alpha1*R1+alpha2*R2)/p; %Like this, I'm not using P for now
Qx = x1-x2;
Qy = y1-y2;
Qz = z1-z2;
KAB = exp(-q*(Qx^2+Qy^2+Qz^2));



end
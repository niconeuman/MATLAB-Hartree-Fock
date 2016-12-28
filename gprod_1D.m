function [E,p,q,P,Q] = gprod_1D(x1,alpha1,x2,alpha2)
p = alpha1+alpha2;
q = (alpha1*alpha2)/p;
P = (alpha1*x1+alpha2*x2)/p; %Like this, I'm not using P for now
Q = x1-x2;
KAB = exp(-q*Q^2);
E = KAB;



end
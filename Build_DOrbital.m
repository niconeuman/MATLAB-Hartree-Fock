function orbital = Build_DOrbital(Ax,Ay,Az,alpha)
%mx = 1, or my = 1 or mz = 1
tmp.x0 = Ax;
tmp.y0 = Ay;
tmp.z0 = Az;
tmp.alpha = alpha;
%L = 2;
%For d functions I don't have the same normalization constant for
%x^2,y^2,z^2, than for xy, xz, yz, and I need to generate a tmp.N vector

prefactor = 2^(2)*(2/pi)^.75*(alpha^(0.75+0.5*(2)));
%xx
%xy
%xz
%yy
%yz
%zz
indices = [2,0,0;1,1,0;1,0,1;0,2,0;0,1,1;0,0,2];
double_factorial_product = zeros(6,1);
for k = 1:6
    double_factorial_product(k) = 1/double_factorial(2*indices(k,1)-1,2*indices(k,2)-1,2*indices(k,3)-1)^0.5;
end
tmp.N = prefactor./double_factorial_product;
tmp.L = 2;
orbital = tmp;


end
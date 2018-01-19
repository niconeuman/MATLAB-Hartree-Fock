function [output,x,index,x_index,Dx] = Interpolated_Boys(n,x,Boys_Table)

%This code will be used once to generate a Table of Boys functions up to
%order 20, x between 0 and 1, and 100 points. That table will be written
%down in this function, and then timings will be compared with the naive
%Boys function

%{
Nlast = 100; %I don't know which is the maximum number possible for x. It
%will depend on g1.alpha, g2.alpha and RPC^2; so probably it can be
%precalculated for a certain system.
Npoints = 1000;
x = linspace(Nlast/Npoints,Nlast,Npoints)';
n = (0:1:30);

[n,x] = meshgrid(n,x);
Boys_Table = gamma(n+0.5).*gammainc(x,n+0.5)./(2*x.^(n+0.5));

n = (0:1:30);Boys_Table = [1./(2*n+1); Boys_Table];
%}
% Nlast = 100;
% Npoints = 200;


%I must define the following variable as global outside the
%Interpolated_Boys function and then call it as global in this function and
%use it.
%The interpolation formula from Helgaker's book is the following:

if (x > 40)
    index = NaN;
    x_index = NaN;
    Dx = NaN;
%     sqrt_pi = 1.772453850905516;
%     %Eq. 9.8.9 of Helgaker's book
%     %F(n,x)= (2n-1)!!/2^(n+1)*sqrt(pi/x^(2n+1));
%     %From wikipedia
%     %(2*k-1)!! = (2k)!/2^k*k! (for k = 1,2,etc)
%     
%     %note: prod(1:n) is at least 3 times faster than factorial(n)
%     if n > 0
%     %output = (prod(1:2*n)/2^n/prod(1:n))/2^(n+1)*sqrt(pi/x^(2*n+1));
%     %output = (prod(1:2*n)/prod(1:n))/2^(2*n+1)*sqrt(pi)/x^(n+0.5);
%     output = (prod(n+1:2*n))/2^(2*n+1)*sqrt_pi/x^(n+0.5);
%     else
%     output = 1/2*sqrt_pi/sqrt(x);
%     end
    
    Prefactor = [0.886226925452758
    443.113462726379e-003
    664.670194089569e-003
    1.66167548522392e+000
    5.81586419828372e+000
    26.1713888922768e+000
    143.942638907522e+000
    935.627152898894e+000
    7.01720364674171e+003
    59.6462309973045e+003
    566.639194474393e+003
    5.94971154198112e+006
    68.4216827327829e+006
    855.271034159787e+006
    11.5461589611571e+009
    167.419304936778e+009
    2.59499922652006e+012
    42.8174872375810e+012
    749.306026657668e+012
    13.8621614931669e+015
    270.312149116754e+015];


%Boys_1_N = Prefactor(N)./x.^(N+0.5);
% for N = n:-1:1
% 
% output(N+1) = Prefactor(N)/x^(N+0.5);
% end
% output(1) = 1/2*sqrt_pi/sqrt(x);



output = Prefactor(n+1)/x^(n+0.5);

    
    
else
xstep = 0.1; %Nlast = 100 /Npoints = 1000       
index = floor(x/xstep)+1; %400 is the number of points, index will give me the previous x_index point, nearest my input x variable
x_index = (index-1)*xstep;
Dx = (x-x_index); %Difference which enters the Taylor expansion
 Dx2 = Dx*Dx;
 Dx3 = Dx2*Dx;
 Dx4 = Dx2*Dx2;
 Dx5 = Dx3*Dx2;

m = n + 1; %column 1 corresponds to n = 0, etc;

output = Boys_Table(index,m)-Boys_Table(index,m+1)*Dx+0.5*Boys_Table(index,m+2)*Dx2-1/6*Boys_Table(index,m+3)*Dx3+1/24*Boys_Table(index,m+4)*Dx4-1/120*Boys_Table(index,m+5)*Dx5;
end
end
function output = Interpolated_Boys_0_small_x(x,Boys_Table)

%Improved version of Interpolated_Boys specific to [ss|ss] type integrals

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
%Nlast = 100;
%Npoints = 200;


%I must define the following variable as global outside the
%Interpolated_Boys function and then call it as global in this function and
%use it.
%The interpolation formula from Helgaker's book is the following:

%if (x > 25)
    %old sqrt_pi = 1.772453850905516;
    %sqrt_pi/2 = 0.886226925452758;
    %Eq. 9.8.9 of Helgaker's book
    %F(n,x)= (2n-1)!!/2^(n+1)*sqrt(pi/x^(2n+1));
    %From wikipedia
    %(2*k-1)!! = (2k)!/2^k*k! (for k = 1,2,etc)
%Because for [ss|ss] n = 0 always, this if is not necessary.    
    %note: prod(1:n) is at least 3 times faster than factorial(n)
%     if n > 0
%     %output = (prod(1:2*n)/2^n/prod(1:n))/2^(n+1)*sqrt(pi/x^(2*n+1));
%     %output = (prod(1:2*n)/prod(1:n))/2^(2*n+1)*sqrt(pi)/x^(n+0.5);
%     output = (prod(n+1:2*n))/2^(2*n+1)*sqrt_pi/x^(n+0.5);
%     else
%     output = 1/2*1.772453850905516/sqrt(x);
%     end
%    output = 0.886226925452758/sqrt(x);
%else
xstep = 0.1;       
index = floor(x/xstep)+1;
x_index = (index-1)*xstep;
Dx = (x-x_index); %Difference which enters the Taylor expansion
 Dx2 = Dx*Dx;
 Dx3 = Dx2*Dx;
 Dx4 = Dx2*Dx2;
 Dx5 = Dx3*Dx2;
%  Dx6 = Dx3*Dx3;
%  Dx7 = Dx4*Dx3;
 half_Dx2 = .5*Dx2;
 sixth_Dx3 = .166666666666667*Dx3;
 twenty_fourth_Dx4 = 4.166666666666666e-2*Dx4;
 hundred_twentieth_Dx5 = 8.333333333333334e-3*Dx5;
%  seven_hundred_twentieth_Dx6 = 1.388888888888889e-3*Dx6;
%  Dx7_over_5040 = 1.984126984126984e-04*Dx7;
 
%m = n + 1; %column 1 corresponds to n = 0, etc;
%n = 0 always!

output = Boys_Table(index,1)-Boys_Table(index,2)*Dx+Boys_Table(index,3)*half_Dx2-Boys_Table(index,4)*sixth_Dx3+Boys_Table(index,5)*twenty_fourth_Dx4-Boys_Table(index,6)*hundred_twentieth_Dx5;%...
        %+Boys_Table(index,6)*seven_hundred_twentieth_Dx6-Boys_Table(index,7)*Dx7_over_5040;
%end
end
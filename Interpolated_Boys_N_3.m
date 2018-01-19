function output = Interpolated_Boys_N_3(n,x,Boys_Table)
% 
%  assert(isa(x,'double'));
%  assert(isa(n,'int32'));
%  assert(all(size(Boys_Table)== [4001 11]));
%  assert(isa(Boys_Table,'double'));
% x = fi(x);
% Boys_Table = fi(Boys_Table);
% n = fi(n);
%This function is similar to Interpolated_Boys, but it calculates at the
%same time, in a vectorized fashion, several Boys functions of order 0 to N

%    Nlast = 100;
%    Npoints = 200;
%N is a vector that goes from 1 to n.


%I must define the following variable as global outside the
%Interpolated_Boys function and then call it as global in this function and
%use it.
%The interpolation formula from Helgaker's book is the following:

if (x > 60)
    
    %sqrt_pi = 1.772453850905516;
    
%The (prod(n+1:2*n)) part generates a vector and multiplies
%its components. So to use it on a vector, which should generate a matrix,
%seems to require for loops. So it is better to precompute it and store it
%on a vector.
%I precompute the quantity (prod(n+1:2*n))/2^(2*n+1)*sqrt_pi
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

%output = zeros(n+1,1);
for N = n+1:-1:1

output(N) = Prefactor(N)/x^(n+0.5);
end


else
%    Nlast = 400;
%    Npoints = 4000;
xstep = 0.1; %Nlast/Npoints
    
index = floor(x/xstep)+1; %400 is the number of points, index will give me the previous x_index point, nearest my input x variable
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
% M = (1:n+1);
% 
% Boys_row = Boys_Table(index,:);
% output = Boys_row(M)-Boys_row(M+1)*Dx+(0.5*Dx2)*Boys_row(M+2)...
%         -(Dx3/6)*Boys_row(M+3)+(Dx4/24)*Boys_row(M+4)...
%         -(Dx5/120)*Boys_row(M+5);
% %output = output';

%This for cycle, with initialization, is 10 times faster than the vector
%code above. Keep in mind!!!!!!!!!!!!!!!!

 output = zeros(n+1,1);
for m = 1:n+1
 output(m) = Boys_Table(index,m)-Boys_Table(index,m+1)*Dx+Boys_Table(index,m+2)*half_Dx2-Boys_Table(index,m+3)*sixth_Dx3+Boys_Table(index,m+4)*twenty_fourth_Dx4-Boys_Table(index,m+5)*hundred_twentieth_Dx5;%...
        %+Boys_Table(index,m+6)*seven_hundred_twentieth_Dx6-Boys_Table(index,m+7)*Dx7_over_5040;
end
end


end
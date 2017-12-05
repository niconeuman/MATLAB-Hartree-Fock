function [output,x,index,x_index,Dx] = Interpolated_Boys_N(n,x,Boys_Table)
%This function is similar to Interpolated_Boys, but it calculates at the
%same time, in a vectorized fashion, several Boys functions of order 0 to N
    Nlast = 100;
    Npoints = 200;

%N is a vector that goes from 1 to n.
N = (1:n)';

%I must define the following variable as global outside the
%Interpolated_Boys function and then call it as global in this function and
%use it.
%The interpolation formula from Helgaker's book is the following:

if (x > Nlast)
    index = NaN;
    x_index = NaN;
    Dx = NaN;
    sqrt_pi = 1.772453850905516;
    
Boys_zero = 1/2*sqrt_pi/sqrt(x);
%The (prod(n+1:2*n)) part generates a vector and multiplies
%its components. So to use it on a vector, which should generate a matrix,
%seems to require for loops. So it is better to precompute it and store it
%on a vector.
%I precompute the quantity (prod(n+1:2*n))/2^(2*n+1)*sqrt_pi
Prefactor = [443.113462726379e-003
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


Boys_1_N = Prefactor(N)./x.^(N+0.5);
output = [Boys_zero; Boys_1_N];
    
else

    
index = floor(x*Npoints/Nlast)+1; %400 is the number of points, index will give me the previous x_index point, nearest my input x variable
x_index = (index-1)*Nlast/Npoints;
Dx = (x-x_index); %Difference which enters the Taylor expansion
 Dx2 = Dx*Dx;
 Dx3 = Dx2*Dx;
 Dx4 = Dx2*Dx2;
 Dx5 = Dx3*Dx2;

%m = n + 1; %column 1 corresponds to n = 0, etc;
M = (1:n+1);

output = Boys_Table(index,M)-Boys_Table(index,M+ones(size(M)))*Dx+0.5*Boys_Table(index,M+2*ones(size(M)))*Dx2...
        -(Dx3/6)*Boys_Table(index,M+3*ones(size(M)))+1/24*Boys_Table(index,M+4*ones(size(M)))*Dx4...
        -(Dx5/120)*Boys_Table(index,M+5*ones(size(M)));
output = output';
end


end
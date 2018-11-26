function output = vrrBoys_Table_3(x,Boys_Table)
%Jan 4th 2018
%This function calculates the interpolated Boys function up to order 2,
%using the Boys_Table (currently up to x = 400)

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


output = [Boys_Table(index,1)-Boys_Table(index,2)*Dx+Boys_Table(index,3)*half_Dx2-Boys_Table(index,4)*sixth_Dx3+Boys_Table(index,5)*twenty_fourth_Dx4-Boys_Table(index,6)*hundred_twentieth_Dx5
           Boys_Table(index,2)-Boys_Table(index,3)*Dx+Boys_Table(index,4)*half_Dx2-Boys_Table(index,5)*sixth_Dx3+Boys_Table(index,6)*twenty_fourth_Dx4-Boys_Table(index,7)*hundred_twentieth_Dx5
           Boys_Table(index,3)-Boys_Table(index,4)*Dx+Boys_Table(index,5)*half_Dx2-Boys_Table(index,6)*sixth_Dx3+Boys_Table(index,7)*twenty_fourth_Dx4-Boys_Table(index,8)*hundred_twentieth_Dx5
           Boys_Table(index,4)-Boys_Table(index,5)*Dx+Boys_Table(index,6)*half_Dx2-Boys_Table(index,7)*sixth_Dx3+Boys_Table(index,8)*twenty_fourth_Dx4-Boys_Table(index,9)*hundred_twentieth_Dx5];
% output = [Boys_Table(index,int8(1))-Boys_Table(index,int8(2))*Dx+Boys_Table(index,int8(3))*half_Dx2-Boys_Table(index,int8(4))*sixth_Dx3+Boys_Table(index,int8(5))*twenty_fourth_Dx4-Boys_Table(index,int8(6))*hundred_twentieth_Dx5
%           Boys_Table(index,int8(2))-Boys_Table(index,int8(3))*Dx+Boys_Table(index,int8(4))*half_Dx2-Boys_Table(index,int8(5))*sixth_Dx3+Boys_Table(index,int8(6))*twenty_fourth_Dx4-Boys_Table(index,int8(7))*hundred_twentieth_Dx5
%           Boys_Table(index,int8(3))-Boys_Table(index,int8(4))*Dx+Boys_Table(index,int8(5))*half_Dx2-Boys_Table(index,int8(6))*sixth_Dx3+Boys_Table(index,int8(7))*twenty_fourth_Dx4-Boys_Table(index,int8(8))*hundred_twentieth_Dx5];
% output = double(output);

end
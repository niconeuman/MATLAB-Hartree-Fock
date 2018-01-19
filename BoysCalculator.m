function Boys_Table = BoysCalculator(m,xmax,xstep,dt)

%This function calculates numerically a Boys_Table from order 0 to m,
%xvalues from 0 to xmax, and xstep separation.

Boys_Table = zeros(xmax/xstep+1,m+1);

%It uses the exact integral expression
%Fm(x) = Int(0 to 1) t^(2*m)*exp(-x*t^2)*dt;


t = linspace(dt/2,1-dt/2,1/dt);

for row = 1:size(Boys_Table,1)
    for col = 1:size(Boys_Table,2)
        
        x = (row-1)*xstep;
        m = (col-1);
        
        Integrand = (t.^(2*m).*exp(-x*t.^2))*dt;
        Boys_x_m = sum(Integrand);
        
        Boys_Table(row,col) = Boys_x_m;
        
    end
    
end


end
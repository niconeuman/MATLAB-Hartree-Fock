function output = Boys(n,x)
if x == 0
    output = 1/(2*n+1);
else
    output = gamma(n+0.5)*gammainc(x,n+0.5)/(2*x^(n+0.5));
end


end
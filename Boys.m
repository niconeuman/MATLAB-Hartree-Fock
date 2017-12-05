function output = Boys(n,x)
%Boys function is costly to evaluate.
%Suggestion in Helgaker's book (Eq. 9.8.12)
%is to preevaluate Boys functions of order n for several points in the x
%interval from 0 to 1.
%Then interpolate for a desired x value using a polynomial formula.

if x == 0
    output = 1/(2*n+1);
else
    output = gamma(n+0.5)*gammainc(x,n+0.5)/(2*x^(n+0.5));
end


end
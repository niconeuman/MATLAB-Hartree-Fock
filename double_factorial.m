function df = double_factorial(mx,my,mz)
if mx == -1
    dfx = 1;
elseif mx == 0
    dfx = 1;
elseif mx == 1
    dfx = 1;
elseif mx == 2
    dfx = 2;
end

if my == -1
    dfy = 1;
elseif my == 0
    dfy = 1;
elseif my == 1
    dfy = 1;
elseif my == 2
    dfy = 2;
end

if mz == -1
    dfz = 1;
elseif mz == 0
    dfz = 1;
elseif mz == 1
    dfz = 1;
elseif mz == 2
    dfz = 2;
end
df = dfx*dfy*dfz;
end
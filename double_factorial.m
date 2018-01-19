function df = double_factorial(mx,my,mz)
if mx == -1
    dfx = 1;
elseif mx == 0
    dfx = 1;
elseif mx == 1
    dfx = 1;
elseif mx == 3
    dfx = 3;
elseif mx == 5
    dfx = 15;
elseif mx == 7
    dfx = 105;
end

if my == -1
    dfy = 1;
elseif my == 0
    dfy = 1;
elseif my == 1
    dfy = 1;
elseif my == 3
    dfy = 3;
elseif my == 5
    dfy = 15;
elseif my == 7
    dfy = 105;
end

if mz == -1
    dfz = 1;
elseif mz == 0
    dfz = 1;
elseif mz == 1
    dfz = 1;
elseif mz == 3
    dfz = 3;
elseif mz == 5
    dfz = 15;
elseif mz == 7
    dfz = 105;
end
df = dfx*dfy*dfz;
end
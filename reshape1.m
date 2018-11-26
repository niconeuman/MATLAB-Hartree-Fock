function OutMat = int_reshape(InMat,L1in,L2in,L1out,L2out,unique)

%This function takes a matrix of dimension (DimL1in x DimL2in) and reshapes
%it into a matrix of dimensions (DimL1out x DimL2out), where L1out =
%L1in+/-1 and L2out = L2in-/+1

%It currently will work only in the first two dimensions

DimL1in = (L1in+1)*(L1in+2)/2;
DimL2in = (L2in+1)*(L2in+2)/2;
DimL1out = (L1out+1)*(L1out+2)/2;
DimL2out = (L2out+1)*(L2out+2)/2;

if (DimL2in == DimL2out)

    OutMat = [InMat(1:DimL1in,:) InMat(unique(1:DimL1in),(end-L2in+1:end)) InMat(unique(1:DimL1in)+1,end)];

elseif (DimL1in == DimL1out)
    %code missing
    OutMat = InMat;
end

end
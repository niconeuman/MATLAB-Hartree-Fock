function OutMat = int_expand_old(InMat,Vec,L,Dim)

%Vec should be a 3x1 column vector
%InMat is a Dim1xDim2 column rectangular matrix
%L is the angular momentum corresponding to the dimension of InMat

if (Dim == 1)
    OutMat = [Vec(1)*InMat;Vec(2)*InMat(end-L:end,:,:,:);Vec(3)*InMat(end,:,:,:)];
elseif (Dim == 2)
    if L == 5
        OutMat = [Vec(1)*InMat(:,1,:,:) Vec(2)*InMat(:,1,:,:) Vec(3)*InMat(:,1,:,:) Vec(2)*InMat(:,2,:,:) Vec(3)*InMat(:,2,:,:) Vec(3)*InMat(:,3,:,:)];           
    else
        OutMat = [Vec(1)*InMat Vec(2)*InMat(:,end-L:end,:,:) Vec(3)*InMat(:,end,:,:)];
    end
end

end
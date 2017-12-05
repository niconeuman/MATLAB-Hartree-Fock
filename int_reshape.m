function [OutMat,DimL1in,DimL2in,DimL1out,DimL2out] = int_reshape(InMat,L1in,L2in,Dim,unique)

%This function takes a matrix of dimension (DimL1in x DimL2in) and reshapes
%it into a matrix of dimensions (DimL1out x DimL2out), where L1out =
%L1in+/-1 and L2out = L2in-/+1

%It currently will work only in the first two dimensions

DimL1in = (L1in+1)*(L1in+2)/2;
DimL2in = (L2in+1)*(L2in+2)/2;

% if L1in == 3 && L2in == 1
%     L1out = L1in-1;
%     L2out = L2in+1;
%     DimL1out = (L1out+1)*(L1out+2)/2;
%     DimL2out = (L2out+1)*(L2out+2)/2;
%     OutMat = [InMat(1:3,1,:,:) InMat(unique(1:3),1,:,:) InMat(unique(1:3)+1,1,:,:)  InMat(unique(1:3),2,:,:)    InMat(unique(1:3),3,:,:)    InMat(unique(1:3)+1,3,:,:)
%               InMat(2,2,:,:)   InMat(4,2,:,:)           InMat(unique(4)+1,1,:,:)    InMat(unique(4),2,:,:)    InMat(unique(4),3,:,:)    InMat(unique(4)+1,3,:,:)
%               InMat(2,3,:,:)   InMat(4,3,:,:)           InMat(unique(5)+1,1,:,:)    InMat(unique(5),2,:,:)    InMat(unique(5),3,:,:)    InMat(unique(5)+1,3,:,:)
%               InMat(3,3,:,:)   InMat(5,3,:,:)           InMat(unique(6)+1,1,:,:)    InMat(unique(6),2,:,:)    InMat(unique(6),3,:,:)    InMat(unique(6)+1,3,:,:)];
% %       OutMat(2,2,:,:) = 0;
% %       OutMat(3,3,:,:) = 0;
% %       OutMat(5,5,:,:) = 0;
% elseif L1in == 4 && L2in == 0
%     OutMat = [InMat(1:10,:,:,:) InMat(unique(1:10),:,:,:) InMat(unique(1:10)+1,:,:,:)];
% %       OutMat(4,1,:,:) = 0;
% %       OutMat(5,1,:,:) = 0;
% %       OutMat(6,1,:,:) = 0;
% %       OutMat(8,1,:,:) = 0;
% %       OutMat(9,1,:,:) = 0;
% else

    if (Dim == 1)
    L1out = L1in-1;
    L2out = L2in+1;
    DimL1out = (L1out+1)*(L1out+2)/2;
    DimL2out = (L2out+1)*(L2out+2)/2;
    OutMat = [InMat(1:DimL1out,:,:,:) InMat(unique(1:DimL1out),(end-L2out+1:end),:,:) InMat(unique(1:DimL1out)+1,end,:,:)];

    elseif (Dim == 2)
    L1out = L1in+1;
    L2out = L2in-1;
    DimL1out = (L1out+1)*(L1out+2)/2;
    DimL2out = (L2out+1)*(L2out+2)/2;
    OutMat = [InMat(:,1:DimL2out);InMat((end-L1out+1:end),unique(1:DimL2out),:,:);InMat(end,unique(1:DimL2out)+1,:,:)];
    end
% end

end
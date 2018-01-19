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
%       OutMat(2,2,:,:) = 0;
%       OutMat(3,3,:,:) = 0;
%       OutMat(5,5,:,:) = 0;
if L1in == 2 && L2in == 0
    if (size(InMat,2) == 1 && size(InMat,3) == 1 && size(InMat,4) == 1)
     OutMat = [InMat(1) InMat(2) InMat(3)
               InMat(2) InMat(4) InMat(5)
               InMat(3) InMat(5) InMat(6)];
    else
     OutMat = [InMat(1,1,:,:) InMat(2,1,:,:) InMat(3,1,:,:)
               InMat(2,1,:,:) InMat(4,1,:,:) InMat(5,1,:,:)
               InMat(3,1,:,:) InMat(5,1,:,:) InMat(6,1,:,:)];
    end
%      OutMat = [InMat(1,1,:,:) InMat(2,1,:,:) InMat(3,1,:,:)
%                InMat(2,1,:,:) InMat(4,1,:,:) InMat(5,1,:,:)
%                InMat(3,1,:,:) InMat(5,1,:,:) InMat(6,1,:,:)];
    
     %OutMat = permute(OutMat,[2 1 3 4]);      
elseif L1in == 3 && L2in == 0
        if (size(InMat,2) == 1 && size(InMat,3) == 1 && size(InMat,4) == 1)
     OutMat = [InMat(1) InMat(2) InMat(3)
               InMat(2) InMat(4) InMat(5)
               InMat(3) InMat(5) InMat(6)
               InMat(4) InMat(7) InMat(8)
               InMat(5) InMat(8) InMat(9)
               InMat(6) InMat(9) InMat(10)];
        else
     OutMat = [InMat(1,1,:,:) InMat(2,1,:,:) InMat(3,1,:,:)
               InMat(2,1,:,:) InMat(4,1,:,:) InMat(5,1,:,:)
               InMat(3,1,:,:) InMat(5,1,:,:) InMat(6,1,:,:)
               InMat(4,1,:,:) InMat(7,1,:,:) InMat(8,1,:,:)
               InMat(5,1,:,:) InMat(8,1,:,:) InMat(9,1,:,:)
               InMat(6,1,:,:) InMat(9,1,:,:) InMat(10,1,:,:)]; 
        end
%      OutMat = [InMat(1,1,:,:) InMat(2,1,:,:) InMat(3,1,:,:)
%                InMat(2,1,:,:) InMat(4,1,:,:) InMat(5,1,:,:)
%                InMat(3,1,:,:) InMat(5,1,:,:) InMat(6,1,:,:)
%                InMat(4,1,:,:) InMat(7,1,:,:) InMat(8,1,:,:)
%                InMat(5,1,:,:) InMat(8,1,:,:) InMat(9,1,:,:)
%                InMat(6,1,:,:) InMat(9,1,:,:) InMat(10,1,:,:)];    
          
elseif L1in == 4 && L2in == 0
    OutMat = [InMat(1:10,:,:,:) InMat(unique(1:10),:,:,:) InMat(unique(1:10)+1,:,:,:)];
elseif L1in == 3 && L2in == 1
%      OutMat = [InMat(1,1,:,:) InMat(1,2,:,:) InMat(1,3,:,:) InMat(2,2,:,:) InMat(2,3,:,:) InMat(3,3,:,:)
%                InMat(2,1,:,:) InMat(2,2,:,:) InMat(2,3,:,:) InMat(4,2,:,:) InMat(4,3,:,:) InMat(5,3,:,:)
%                InMat(3,1,:,:) InMat(3,2,:,:) InMat(3,3,:,:) InMat(5,2,:,:) InMat(5,3,:,:) InMat(6,3,:,:)
%                InMat(4,1,:,:) InMat(4,2,:,:) InMat(4,3,:,:) InMat(7,2,:,:) InMat(7,3,:,:) InMat(8,3,:,:)
%                InMat(5,1,:,:) InMat(5,2,:,:) InMat(5,3,:,:) InMat(8,2,:,:) InMat(8,3,:,:) InMat(9,3,:,:)
%                InMat(6,1,:,:) InMat(6,2,:,:) InMat(6,3,:,:) InMat(9,2,:,:) InMat(9,3,:,:) InMat(10,3,:,:)];
% %     %If I want [dxx dxy] (1,2), I need [fxxy px] (2,1)
    OutMat = [InMat(1,1,:,:) InMat(2,1,:,:) InMat(3,1,:,:) InMat(2,2,:,:) InMat(3,2,:,:) InMat(3,3,:,:)
              InMat(2,1,:,:) InMat(4,1,:,:) InMat(5,1,:,:) InMat(4,2,:,:) InMat(5,2,:,:) InMat(5,3,:,:)
              InMat(3,1,:,:) InMat(5,1,:,:) InMat(6,1,:,:) InMat(5,2,:,:) InMat(6,2,:,:) InMat(6,3,:,:)
              InMat(4,1,:,:) InMat(7,1,:,:) InMat(8,1,:,:) InMat(7,2,:,:) InMat(8,2,:,:) InMat(8,3,:,:)
              InMat(5,1,:,:) InMat(8,1,:,:) InMat(9,1,:,:) InMat(8,2,:,:) InMat(9,2,:,:) InMat(9,3,:,:)
              InMat(6,1,:,:) InMat(9,1,:,:) InMat(10,1,:,:) InMat(9,2,:,:) InMat(10,2,:,:) InMat(10,3,:,:)];          
%      OutMat = [InMat(1:3,1,:,:) InMat(unique(1:3),1,:,:) InMat(unique(1:3)+1,1,:,:)  InMat(unique(1:3),2,:,:)    InMat(unique(1:3),3,:,:)    InMat(unique(1:3)+1,3,:,:)
%                InMat(2,2,:,:)   InMat(4,2,:,:)           InMat(unique(4)+1,1,:,:)    InMat(unique(4),2,:,:)    InMat(unique(4),3,:,:)    InMat(unique(4)+1,3,:,:)
%                InMat(2,3,:,:)   InMat(4,3,:,:)           InMat(unique(5)+1,1,:,:)    InMat(unique(5),2,:,:)    InMat(unique(5),3,:,:)    InMat(unique(5)+1,3,:,:)
%                InMat(3,3,:,:)   InMat(5,3,:,:)           InMat(unique(6)+1,1,:,:)    InMat(unique(6),2,:,:)    InMat(unique(6),3,:,:)    InMat(unique(6)+1,3,:,:)];
              
end          
%       OutMat(4,1,:,:) = 0;
%       OutMat(5,1,:,:) = 0;
%       OutMat(6,1,:,:) = 0;
%       OutMat(8,1,:,:) = 0;
%       OutMat(9,1,:,:) = 0;
% else
% if L1in == 3 && L2in == 1
%     OutMat = [InMat(1:6,1,:,:) InMat(1:6,2,:,:) InMat(1:6,3,:,:) InMat([2,4,5,7,8,9],2,:,:) InMat([2,4,5,7,8,9],3,:,:) InMat([3,5,6,8,9,10],3,:,:)];
% elseif L1in == 2 && L2in == 0
%     OutMat = [InMat(1,1,:,:) InMat(2,1,:,:) InMat(3,1,:,:)
%               InMat(2,1,:,:) InMat(4,1,:,:) InMat(5,1,:,:)
%               InMat(3,1,:,:) InMat(5,1,:,:) InMat(6,1,:,:)];
        
% else
%     if (Dim == 1)
%     L1out = L1in-1;
%     L2out = L2in+1;
%     DimL1out = (L1out+1)*(L1out+2)/2;
%     DimL2out = (L2out+1)*(L2out+2)/2;
%     OutMat = [InMat(1:DimL1out,:,:,:) InMat(unique(1:DimL1out),(end-L2out+1:end),:,:) InMat(unique(1:DimL1out)+1,end,:,:)];
% 
%     elseif (Dim == 2)
%     L1out = L1in+1;
%     L2out = L2in-1;
%     DimL1out = (L1out+1)*(L1out+2)/2;
%     DimL2out = (L2out+1)*(L2out+2)/2;
%     OutMat = [InMat(:,1:DimL2out);InMat((end-L1out+1:end),unique(1:DimL2out),:,:);InMat(end,unique(1:DimL2out)+1,:,:)];
%     end
%  end

end
function gabcd = ascs(a,b,c,d,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1,L2,L3,L4,Lmax,order,Boys_Table,nz,unique)
%[as|cs] = -p/q*[(a+1)s|(c-1)s]+(RQC+p/q*RPA)*[as|(c-1)s]
%          +a_i/2/q*[(a-1)s|(c-1)s]+c_i/2/q*[as|(c-2)s]

%For example if [as|cs] is [ds|ps] (6 by 3), then the first term will be
%[fs|ss] (10 x 1 reshaped as 18 x 1), the second term expand((RQC+p/q*RPA),[ds|ss])
%the 3rd term will be TwiceExp([ps|ss]).*nz{3} (This will give a [fs|ss]
%type matrix which will have to be reshaped as (18 x 1)
%The 4th term will be 0.

if L3 == 1
    
    %[ps|ps] = [3,1,3,1]
    
    %The gLaSLcm1s output for this case will be a (Dim(L1),1,1,1) matrix
    %The gLam1SLcm1S output will be a (Dim(L1-1),1,1,1) matrix
    %[[ps|ss],~,[ss|ss],~] = [3,1,1,1], ,[1,1,1,1]
    [gLaSLcm1s,~,gLam1SLcm1S,~] = vrr_1(a,b,c,d,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1,L2,L3-1,L4,Lmax,order,Boys_Table,nz);
    %[[ds|ss],~,~[ps|ss],~] = [6,1,1,1]
    [gLap1SLcm1s,~,~,~] = vrr_1(a,b,c,d,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1+1,L2,L3-1,L4,Lmax,order,Boys_Table,nz);
    
    %As an example
%     gLap1SLcm1s = [ds|ss] = [xx        ? [xx xy xz
%                             xy            xy yy yz
%                             xz            xz yz zz]
%                             yy
%                             yz
%                             zz]
%   [ds(1:3) ds([2 4 5] ds([3 5 6])];
%La = L1 + 1; %I don't want to change L1 
DimL1 = (L1+1)*(L1+2)/2; %L1 = 1, DimL1 = 3
%Dima = (La+1)*(La+2)/2; %La = 2 in this case, %Dima = 6
Term1 = [gLap1SLcm1s(1:DimL1) gLap1SLcm1s(unique(1:DimL1)) gLap1SLcm1s(unique(1:DimL1)+1)];
    %Another example
    %fs = [xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz]' -> 
    %[fs(1:6) fs([2 4 5 7 8 9]) fs([3 5 6 8 9 10])]
   
    
    %[xxX xxY xxZ                
    % xyX                        
    % xzX         =          
    % yyX
    % yzX
    % zzX    zzZ]
    
    %L1 = 3, DimL1 = 6
    %V_fs = [V_L1_p1_L2_m1(1:6),V_L1_p1_L2_m1([2 4 5 7 8 9]),V_L1_p1_L2_m1([3 5 6 8 9 10])];
    %V_fs = [gLap1SLcm1s(1:DimL1) gLap1SLcm1s([DimL1-1,DimL1+1:2*DimL1-1]) gLap1SLcm1s([DimL1,DimL1+2:2*DimL1])];

%Finally Term1 = [3,3,1,1] %Remember I have permuted dimensions to make
%matrix operations easier.
%Term1 = [gLap1SLcm1s(1:end)];
    
    Vec2 = RQC+p/q*RPA;
    %Expand(Vec2,gLaSLcm1s); % ([3x1],[3x1])
    ExpTerm2 = [Vec2(1)*gLaSLcm1s;Vec2(2)*gLaSLcm1s(end-L1:end);Vec2(3)*gLaSLcm1s(end)];
    %The following line is reshaping the previous 6 x 1 matrix into a 3 x 3
    %matrix
    Term2 = [ExpTerm2(1:DimL1) ExpTerm2(unique(1:DimL1)) ExpTerm2(unique(1:DimL1)+1)];
    
    PrevTerm3 = gLam1SLcm1S/2/q;
    if (L1 == 1)
    %This is because PrevTerm3 = [ss|ss]
    TwiceExpTerm3 = nz{L1+1}.*(PrevTerm3*ones(6,1));
    else %L1 ==2, etc
    ExpTerm3 = [PrevTerm3;PrevTerm3(end-(L1-1)+1:end);PrevTerm3(end)];
    TwiceExpTerm3 = nz{L1+1}.*[ExpTerm3;ExpTerm3(end-L1+1:end);ExpTerm3(end)];
    end
    
    Term3 = [TwiceExpTerm3(1:DimL1) TwiceExpTerm3(unique(1:DimL1)) TwiceExpTerm3(unique(1:DimL1)+1)];
    %Term3 = [PrevTerm3;PrevTerm3(end-L1+1:end);PrevTerm3(end)];
    gLaSLcs = -p/q*Term1+Term2+Term3;%+nz{L3}/2/q*[as|(c-2i)s];
    
             
    gabcd = gLaSLcs;
elseif L3 == 2
    %[as|cs] = -p/q*[(a+1)s|(c-1)s]+(RQC+p/q*RPA)*[as|(c-1)s]
%             +a_i/2/q*[(a-1)s|(c-1)s]+c_i/2/q*[as|(c-2)s]
    
    %L1 >= L3, has been taken care of at Build_ERI_2
    %Cases
    %[ds|ds] = [6,1,6,1]
    
    %[[ds|ss],~,~,~] = [6,1,1,1]
    [gLaSLcm2S,~,~,~] = vrr_1(a,b,c,d,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1,L2,L3-2,L4,Lmax,order,Boys_Table,nz);
    %[[ps|ps],~,~,~] = [3,1,3,1] 
    [gLam1SLcm1S,~,~,~] = vrr_1(a,b,c,d,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1-1,L2,L3-1,L4,Lmax,order,Boys_Table,nz);
    %[[fs|ps],~,~[ds|ps],~] = [10,1,3,1],[6,1,3,1]
    [gLap1SLcm1S,~,gLaSLcm1S,~] = vrr_1(a,b,c,d,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1+1,L2,L3-1,L4,Lmax,order,Boys_Table,nz);

    DimL1 = (L1+1)*(L1+2)/2; %L1 = 2, DimL1 = 6
    
    %The size of this matrix should be [6x6]
    %Taken from HRR_Nuc
    %V_fp = [V_L1_p1_L2_m1(1:6,:),V_L1_p1_L2_m1([2 4 5 7 8 9],2:3),V_L1_p1_L2_m1([3 5 6 8 9 10],3)];
    Term1 = [gLap1SLcm1S(1:DimL1,:) gLap1SLcm1S(unique(1:DimL1),2:3) gLap1SLcm1S(unique(1:DimL1)+1,3)];

    Vec2 = RQC+p/q*RPA;
    %Expand(Vec2,gLaSLcm1s); % ([3x1],[6x3]) -> [6x6] (So I expand in
    %direction 2)
    %Example, also from HRR_Nuc:
    %[RAB(1)*V_dp,RAB(2)*V_dp(:,2:3);RAB(3)*V_dp(:,3)]
    %Term2 = [Vec2(1)*gLaSLcm1S;Vec2(2)*gLaSLcm1S(:,2:3);Vec2(3)*gLaSLcm1S(:,3)];
    Term2 = [Vec2(1)*gLaSLcm1S Vec2(2)*gLaSLcm1S(:,2:3) Vec2(3)*gLaSLcm1S(:,3)];
    
    %Term 3
    %[ps|ps]
    %+a_i/2/q*[(a-1)s|(c-1)s]
    
    %gLam1SLcm1S = [3,3,1,1]
    PrevTerm3 = gLam1SLcm1S/2/q;
    %To get a 6x6 term I should expand once in each dimension
    
    %With the following line, I get a 6x3 matrix
    ExpTerm3 = [PrevTerm3;PrevTerm3(end-(L1-1)+1:end,:);PrevTerm3(end,:)];
    
    %The following should give me a 6 x 6 matrix, but I don't know if it's
    %OK
    
    %TwiceExpTerm3 = nz{L1}*nz{L3}'.*[ExpTerm3;ExpTerm3(:,end-L1+1:end);ExpTerm3(:,end)];
    %Shouldn't it be a row concatenation?
    TwiceExpTerm3 = nz{L1}*nz{L3}'.*[ExpTerm3 ExpTerm3(:,end-L1+1:end) ExpTerm3(:,end)];
    
    %I have to be careful with nz{L1}*nz{L3}'. The size seems correct,
    %because L1 = 2, but I don't know if the logic is OK
    
    %Term 4
    %+c_i/2/q*[as|(c-2)s]
    %[[ds|ss],~,~,~] = [6,1,1,1]
    PrevTerm4 = gLaSLcm2S/2/q;
    %This is because PrevTerm4 = [ds|ss]
    TwiceExpTerm4 = nz{L3}.*(PrevTerm4*ones(1,6));
    
    gLaSLcS = -p/q*Term1+Term2+TwiceExpTerm3+TwiceExpTerm4;%+nz{L3}/2/q*[as|(c-2i)s];
    
    gabcd = gLaSLcS;
elseif L3 == 3
else %L3 > 3
end


end
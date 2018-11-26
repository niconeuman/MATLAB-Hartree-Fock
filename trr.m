function gLaSLcs = trr(a,b,c,d,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1,L2,L3,L4,Lmax,order,Boys_Table,nz)

%April 12th 2017

%This function increases angular momentum in the c function of an ERI
%trr = transfer recursion relation

%The recursion formula is
%[as|(c+1i)s] = -p/q*[(a+1i)s|cs]+((Qi-Ci)+p/q*(Pi-Ai))*[as|cs]
%                +nz{L1}/2/q*[(a-1i)s|cs]+nz{L3}/2/q*[as|(c-1i)s];
%or
%[as|Cs] = -p/q*[(a+1i)s|(c-1i)s]+((Qi-Ci)+p/q*(Pi-Ai))*[as|(c-1i)s]
%                +nz{L1}/2/q*[(a-1i)s|(c-1i)s]+nz{L3}/2/q*[as|(c-2i)s];
%Careful with the nz vectors, I should check that array dimensions are
%always correct.

if L3 == 1
    
    %The gLaSLcm1s output for this case will be a (Dim(L1),1,1,1) matrix
    %The gLam1SLcm1S output will be a (Dim(L1-1),1,1,1) matrix
    [gLaSLcm1s,~,gLam1SLcm1S,~] = vrr_1(a,b,c,d,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1,L2,0,L4,Lmax,order,Boys_Table,nz);
    [gLap1SLcm1s,~,~,~] = vrr_1(a,b,c,d,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1+1,L2,0,L4,Lmax,order,Boys_Table,nz);

    Vec2 = RQC+p/q*RPA;
    Term2 = [Vec2(1)*gLaSLcm1s;Vec2(2)*gLaSLcm1s(end-L3+1:end);Vec2(3)*gLaSLcm1s(end)];
    
    PrevTerm3 = gLam1SLcm1S/2/q;
    Term3 = [PrevTerm3;PrevTerm3(end-L1+1:end);PrevTerm3(end)];
    gLaSLcs = -p/q*gLap1SLcm1s+Term2...
                 +nz{L1}.*Term3;%+nz{L3}/2/q*[as|(c-2i)s];
elseif L3 == 2
elseif L3 == 3
else %L3 > 3
end



















end
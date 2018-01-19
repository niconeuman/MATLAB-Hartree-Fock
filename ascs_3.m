function gabcd = ascs_3(aa,ab,ac,ad,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1,L2,L3,L4,Lmax,order,Boys_Table,nz,unique)
%[as|cs] = -p/q*[(a+1)s|(c-1)s]+(RQC+p/q*RPA)*[as|(c-1)s]
%          +a_i/2/q*[(a-1)s|(c-1)s]+c_i/2/q*[as|(c-2)s]

%For example if [as|cs] is [ds|ps] (6 by 3), then the first term will be
%[fs|ss] (10 x 1 reshaped as 18 x 1), the second term expand((RQC+p/q*RPA),[ds|ss])
%the 3rd term will be TwiceExp([ps|ss]).*nz{3} (This will give a [fs|ss]
%type matrix which will have to be reshaped as (18 x 1)
%The 4th term will be 0.
if L3 == 0
    %disp('Is this used');
    gabcd = vrr_1minimal(Kab,Kcd,RPA,RWP,p,q,alpha,RPQ2,L1,order,Boys_Table);

elseif L3 == 1 %[ps|ps], [ds|ps], [fs|ps], etc
    
    %[ps|ps] = [3,1,3,1]
    
    %The gLaSLcm1s output for this case will be a (Dim(L1),1,1,1) matrix
    %The gLam1SLcm1S output will be a (Dim(L1-1),1,1,1) matrix
    %[[ps|ss],~,[ss|ss],~] = [3,1,1,1], ,[1,1,1,1]
                                      %(a,b,c,d,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1,L2,L3,L4,Lmax,order,Boys_Table,nz)
    [gLaSLcm1s,~,gLam1SLcm1s,~] = vrr_1minimal(Kab,Kcd,RPA,RWP,p,q,alpha,RPQ2,L1,order,Boys_Table);
    %[[ds|ss],~,~[ps|ss],~] = [6,1,1,1]
    [gLap1SLcm1s,~,~,~] = vrr_1minimal(Kab,Kcd,RPA,RWP,p,q,alpha,RPQ2,L1,order,Boys_Table);
    
    %gLap1SLcm1s %[ds|ss], [fs|ss], [gs|ss], etc
    %gLaSLcm1s %[ps|ss], [ds|ss], [fs|ss], etc
    
    term1 = -(p/q)*int_reshape(gLap1SLcm1s,L1+1,L3-1,1,unique); %[ps|ps], [ds|ps], [fs|ps], etc
    term2 = int_expand(gLaSLcm1s,(RQC+p/q*RPA),L3-1,2); %[ps|ps], [ds|ps], [fs|ps], etc
    %Apparently in Lindh,1991, the term RQC+p/q*RPA is actually
    %RQC+alpha/q*RPA (have to check it).
    %term2 = int_expand(gLaSLcm1s,(RQC+alpha/q*RPA),L3-1,2);
    

        
        if L1 == 1 %ssss gLam1SLcm1s
            term3 = 1/2/q*eye(3)*gLam1SLcm1s;
        elseif L1 == 2 %psss gLam1SLcm1s
            term3 = 1/2/q*[2*gLam1SLcm1s(1) 0                   0
                           gLam1SLcm1s(2)   gLam1SLcm1s(1)      0
                           gLam1SLcm1s(3)   0                   gLam1SLcm1s(1)
                           0                2*gLam1SLcm1s(2)    0
                           0                gLam1SLcm1s(3)      gLam1SLcm1s(2)
                           0                0                   2*gLam1SLcm1s(3)];            
         elseif L1 == 3 %dsss gLam1SLcm1s
             term3 = 1/2/q*[3*gLam1SLcm1s(1)    0                   0
                            2*gLam1SLcm1s(2)    1*gLam1SLcm1s(1)    0
                            2*gLam1SLcm1s(3)    0                   1*gLam1SLcm1s(1)
                            1*gLam1SLcm1s(4)    2*gLam1SLcm1s(2)    0
                            1*gLam1SLcm1s(5)    1*gLam1SLcm1s(3)    1*gLam1SLcm1s(2)
                            1*gLam1SLcm1s(6)    0                   2*gLam1SLcm1s(3)
                            0                   3*gLam1SLcm1s(4)    0
                            0                   2*gLam1SLcm1s(5)    1*gLam1SLcm1s(4)
                            0                   1*gLam1SLcm1s(6)    2*gLam1SLcm1s(5)
                            0                   0                   3*gLam1SLcm1s(6)];                
        end

    gabcd = term1 + term2 + term3;


elseif L3 == 2 %[ds|ds], [fs|ds], [gs|ds], etc
    %[as|cs] = -p/q*[(a+1)s|(c-1)s]+(RQC+p/q*RPA)*[as|(c-1)s]
    %        +a_i/2/q*[(a-1)s|(c-1)s]+c_i/2/q*[as|(c-2)s]
    
    %For now this code is inefficient.
    %[ps|ps], [ds|ps], [fs|ps], etc:
    gam1scm1s = ascs_3(aa,ab,ac,ad,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1-1,L2,L3-1,L4,Lmax,order,Boys_Table,nz,unique);
    
    %[ds|ss], [fs|ss], [gs|ss], etc:
    gascm2s = ascs_3(aa,ab,ac,ad,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1,L2,L3-2,L4,Lmax,order,Boys_Table,nz,unique);
    
    %[ds|ps], [fs|ps], [gs|ps], etc:
    gascm1s = ascs_3(aa,ab,ac,ad,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1,L2,L3-1,L4,Lmax,order,Boys_Table,nz,unique);
    
    %[fs|ps], [gs|ps], [hs|ps], etc:
    gap1scm1s = ascs_3(aa,ab,ac,ad,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1+1,L2,L3-1,L4,Lmax,order,Boys_Table,nz,unique);
    
    term1 = -(p/q)*int_reshape(gap1scm1s,L1+1,L3-1,1,unique);
    term2 = int_expand(gascm1s,(RQC+p/q*RPA),L3-1,2);
    %Apparently in Lindh,1991, the term RQC+p/q*RPA is actually
    %RQC+alpha/q*RPA (have to check it).
    %term2 = int_expand(gascm1s,(RQC+alpha/q*RPA),L3-1,2);
    

    nz22_term3 =   [2*gam1scm1s(1,1)    0                   0                   0                   0                   0
                    1*gam1scm1s(2,1)    1*gam1scm1s(1,1)    0                   1*gam1scm1s(1,2)    0                   0
                    1*gam1scm1s(3,1)    0                   1*gam1scm1s(1,1)    0                   1*gam1scm1s(1,2)    1*gam1scm1s(1,3)
                    0                   1*gam1scm1s(2,1)    0                   2*gam1scm1s(2,2)    0                   0
                    0                   1*gam1scm1s(3,1)    1*gam1scm1s(2,1)    1*gam1scm1s(3,2)    1*gam1scm1s(2,2)    1*gam1scm1s(2,3)
                    0                   0                   1*gam1scm1s(3,1)    0                   1*gam1scm1s(3,2)    2*gam1scm1s(3,3)];


    term3 = 1/2/q*nz22_term3;    

%Esto parecía andar (puede que no tuviera efecto, para CH4)
    nz22_term4 = [gascm2s(1) 0 0 gascm2s(1) 0 gascm2s(1)
                  gascm2s(2) 0 0 gascm2s(2) 0 gascm2s(2)
                  gascm2s(3) 0 0 gascm2s(3) 0 gascm2s(3)
                  gascm2s(4) 0 0 gascm2s(4) 0 gascm2s(4)
                  gascm2s(5) 0 0 gascm2s(5) 0 gascm2s(5)
                  gascm2s(6) 0 0 gascm2s(6) 0 gascm2s(6)];           

    term4 = 1/2/q*nz22_term4;
    
    gabcd = term1 + term2 + term3 + term4;


elseif L3 == 3
    disp('This should not happen, L3 = 3 in ascs_3');
    %[as|cs] = -p/q*[(a+1)s|(c-1)s]+(RQC+p/q*RPA)*[as|(c-1)s]
    %        +a_i/2/q*[(a-1)s|(c-1)s]+c_i/2/q*[as|(c-2)s]
    
    %For now this code is inefficient.
    %[ps|ps], [ds|ps], [fs|ps], etc:
    gam1scm1s = ascs_3(aa,ab,ac,ad,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1-1,L2,L3-1,L4,Lmax,order,Boys_Table,nz,unique);
    
    %[ds|ss], [fs|ss], [gs|ss], etc:
    gascm2s = ascs_3(aa,ab,ac,ad,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1,L2,L3-2,L4,Lmax,order,Boys_Table,nz,unique);
    
    %[ds|ps], [fs|ps], [gs|ps], etc:
    gascm1s = ascs_3(aa,ab,ac,ad,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1,L2,L3-1,L4,Lmax,order,Boys_Table,nz,unique);
    
    %[fs|ps], [gs|ps], [hs|ps], etc:
    gap1scm1s = ascs_3(aa,ab,ac,ad,Kab,Kcd,RAB,RCD,RPA,RPB,RQC,RQD,RWP,RWQ,p,q,RPQ,alpha,RPQ2,L1+1,L2,L3-1,L4,Lmax,order,Boys_Table,nz,unique);
    
    term1 = -(p/q)*int_reshape(gap1scm1s,L1+1,L3-1,1,unique);
    term2 = int_expand(gascm1s,(RQC+p/q*RPA),L3-1,2);
    %Apparently in Lindh 1991, its alpha, not p
    %term2 = int_expand(gascm1s,(RQC+alpha/q*RPA),L3-1,2);
    prevterm3 = gam1scm1s; %[ps|ps], [ds|ps], [fs|ps], etc
    %This is ones(1,6), because all these matrices have d on the rows
    expterm3 = (nz{L1}*ones(1,3)).*int_expand(prevterm3,ones(3,1),L1-1,1);
    term3 = int_expand(expterm3,ones(3,1),L3-1,2)/2/q;
    
    %I try this
    
%       expterm3 = int_expand(prevterm3,ones(3,1),L1-1,1);
%       term3 = ((nz{L1})*nz{L3}').*int_expand(expterm3,ones(3,1),L3-1,2)/2/q;
    
    %Dim1 = (L1+1)*(L1+2)/2; 
    prevterm4 = gascm2s; %[ds|ss], [fs|ss], [gs|ss], etc
    expterm4 = int_expand(prevterm4,ones(3,1),L3-2,2); %[ds|ps], [fs|ps], [gs|ps], etc
    term4 = (nz{L1}*nz{L3}').*int_expand(expterm4,ones(3,1),L3-1,2)/2/q; %[ds|ds], [fs|ds], [gs|ds], etc
    
    gabcd = term1 + term2 + term3 + term4;
else %L3 > 3
end


end
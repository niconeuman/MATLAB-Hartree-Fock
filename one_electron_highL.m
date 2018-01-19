function [S,T] = one_electron_highL(g1,g2,L1,L2,RPA,RAB,nz,unique)
%December 30th 2016
%Is it worth to use recursion relations for the one electron integrals?
%If I have f-functions (L = 3), I need the cases
%0,0 -------------------------------1 case
%0,1 1,1 1,0 -----------------------3 cases (4 total)
%2,0 2,1 2,2 1,2 0,2 ---------------5 cases (9 total)
%3,0 3,1 3,2 3,3 2,3 1,3 0,3 -------7 cases (16 total)

%I can certainly program all these cases, but when I add g-functions there
%will be 25 cases. And each case will contain all previous data.

%My choice is that one_electron will only deal with s and p functions
%and one_electron_highL will deal recursively with higher angular momentum
%functions

a = g1.alpha;
b = g2.alpha;
p = a+b;
q = a*b/p;
%P = [a*g1.x0 + b*g2.x0,a*g1.y0 + b*g2.y0,a*g1.z0 + b*g2.z0]/p;
R12sq = (g1.x0-g2.x0)^2+(g1.y0-g2.y0)^2+(g1.z0-g2.z0)^2;

if (L1 == 0 && L2 == 0) %[s|s]

    KAB = exp(-q*(RAB(1)^2+RAB(2)^2+RAB(3)^2));
    
    S = KAB*sqrt(pi/(g1.alpha + g2.alpha))^3;
    T = S*(a*b/p*(3-2*a*b/p*R12sq));

elseif (L1 > 0 && L2 == 0)
    [Sam1b,Tam1b] = one_electron_highL(g1,g2,L1-1,L2,RPA,RAB,nz);
    if L1 == 1    %[p|s]
        S = RPA*Sam1b;
        T = RPA*Tam1b+2*q*S;
    elseif L1 == 2 %[d|s]
        [Sam2b,Tam2b] = one_electron_highL(g1,g2,L1-2,L2,RPA,RAB,nz);
        %S = [RPA(1)*Sam1b;RPA(2)*Sam1b(end-L1+1:end);RPA(3)*Sam1b(end)]+nz{2}/2/p*Sam2b;
        %S = int_expand(RPA,Sam1b,L1-1,1)+nz{2}/2/p*Sam2b;
        %STam2b = (Tam2b-2*q*Sam2b)/2/p;
        %T = [RPA(1)*Tam1b;RPA(2)*Tam1b(end-L1+1:end);RPA(3)*Tam1b(end)]+2*q*S+nz{2}*STam2b;
%26/12/207. I try compatibilizing the previous expressions to what is used
%for two electron integrals
        Sterm1 = int_expand(Sam1b,RPA,L1-1,1);
        Sterm2 = 1/2/p*[Sam2b
                        0
                        0
                        Sam2b
                        0
                        Sam2b];
        S = Sterm1+Sterm2;            
        STam2b = (Tam2b-2*q*Sam2b)/2/p;
        Tterm1 = int_expand(Tam1b,RPA,L1-1,1);
        Tterm2 = 2*q*S;
        Tterm3 = [STam2b
                  0
                  0
                  STam2b
                  0
                  STam2b];
        T = Tterm1+Tterm2+Tterm3;
        
    else %[f|s], [g|s], etc
        [Sam2b,Tam2b] = one_electron_highL(g1,g2,L1-2,L2,RPA,RAB,nz);
        ExpSam2b = [Sam2b;Sam2b(end-L1+2:end);Sam2b(end)];
        Exp2Sam2b = [ExpSam2b;ExpSam2b(end-L1+1:end);ExpSam2b(end)];
        S = [RPA(1)*Sam1b;RPA(2)*Sam1b(end-L1+1:end);RPA(3)*Sam1b(end)]+nz{L1}.*Exp2Sam2b/2/p;
        STam2b = (Tam2b-2*q*Sam2b)/2/p;
        ExpSTam2b = [STam2b;STam2b(end-L1+2:end);STam2b(end)];
        Exp2STam2b = [ExpSTam2b;ExpSTam2b(end-L1+1:end);ExpSTam2b(end)];
        T = [RPA(1)*Tam1b;RPA(2)*Tam1b(end-L1+1:end);RPA(3)*Tam1b(end)]+2*q*S+nz{L1}.*Exp2STam2b;
    end

elseif (L1 > 0 && L2 > 0)
    %dummy code so the program doesn't crash due to singularity of the S
    %matrix
    %Dim1 = (L1+1)*(L1+2)/2;
    %Dim2 = (L2+1)*(L2+2)/2;
    %S = eye(Dim1,Dim2);
    %T = eye(Dim1,Dim2);
    
    %[a|b] = [a|b-1]*RPA + [a+1|b-1]
    
    [Sabm1,Tabm1] = one_electron_highL(g1,g2,L1,L2-1,RPA,RAB,nz,unique);
    [Sap1bm1,Tap1bm1] = one_electron_highL(g1,g2,L1+1,L2-1,RPA,RAB,nz,unique);
    %Case [p|p] needs [p|s] and [d|s]
    ExpSabm1 = int_expand(Sabm1,RAB,L2-1,2);
    %termS1 = int_reshape(ExpSabm1,L1+1,L2-1,1,unique);
    termS2 = int_reshape(Sap1bm1,L1+1,L2-1,1,unique);
        
        S = ExpSabm1+termS2;
        %S = termS1+termS2;
    if L2 == 1
        
        term1 = int_reshape(Tap1bm1,L1+1,L2-1,1,unique);
        term2 = int_expand(Tabm1,RAB,L2-1,2);
        term3 = S;
        term4 = termS2;
        term5 = 0;
        [Sam1bm1,~] = one_electron_highL(g1,g2,L1-1,L2-1,RPA,RAB,nz,unique);
        %Sam1bm1 would be an [s|s] scalar or a [p|s] or higher vector
        %prevterm6 = int_expand(Sam1bm1,ones(3,1),L1-1,1);
        %term6 = 1/2/p*int_reshape(nz{L1+1},L1+1,L2-1,1,unique).*int_expand(prevterm6,ones(3,1),L2-1,2);
        if L1 == 1
             term6 = 1/2/p*[Sam1bm1     0       0
                            0           Sam1bm1 0
                            0           0       Sam1bm1];
        elseif L1 == 2
             term6 = 1/2/p*[2*Sam1bm1(1)       0                   0
                            Sam1bm1(2)         Sam1bm1(1)          0
                            Sam1bm1(3)         0                   0
                            0                  2*Sam1bm1(2)        Sam1bm1(1)
                            0                  Sam1bm1(3)          Sam1bm1(2)
                            0                  0                   2*Sam1bm1(3)];
        end
        T = term1 + term2 +2*q*(term3-term4-term5+term6);
    elseif L2 == 2
        %This whole part is not yet debugged.
        disp('One_electron with L2 == 2 not yet implemented');
    [Sabm2,Tabm2] = one_electron_highL(g1,g2,L1,L2-2,RPA,RAB,nz,unique);
    %Sabm2 would be an [d|s] (6 x 1) matrix or higher N x 1 matrix
        term1 = int_reshape(Tap1bm1,L1+1,L2-1,1,unique);
        term2 = int_expand(Tabm1,RAB,L2-1,2);
        term3 = S;
        term4 = termS2;
        prevterm5 = int_expand((nz{L1}.*Sabm2),ones(3,1),L1-2,2); %This would be an 
        term5 = 1/2/p*int_expand(prevterm5,ones(3,1),L1-1,2);
        [Sam1bm1,Tam1bm1] = one_electron_highL(g1,g2,L1-1,L2-1,RPA,RAB,nz,unique);
        %Sam1bm1 would be an [p|p] matrix or a [d|p] or higher N x 3 matrix
        prevterm6 = int_expand(Sam1bm1,ones(3,1),L1-1,1);
        term6 = 1/2/p*int_reshape(nz{L1+1},L1+1,L2-1,1,unique).*int_expand(prevterm6,ones(3,1),L2-1,2);
        T = term1 + term2 +2*q*(term3-term4-term5+term6);
    else %L2 > 2
        
    end

end



end
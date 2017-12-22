function [VAB,nAtom] = Build_Nuclear_Attraction_2(basis,Shell_Doublets,NShell_Doublets,AL,Z,Boys_Table)
nAtom = size(Z,2);
Ncont = Shell_Doublets(end,2); %This is the final mu_end
VAB = zeros(Ncont,Ncont);
nz = cell(1,10);
nz{2} = [1;0;0;1;0;1];
nz{3} = [2;1;1;0;0;0;2;1;0;2];
nz{4} = [3;2;2;1;1;1;0;0;0;0;3;2;1;0;3];
nz{5} = [4;3;3;2;2;2;1;1;1;1;0;0;0;0;0;4;3;2;1;0;4];
nz{6} = [5;4;4;3;3;3;2;2;2;2;1;1;1;1;1;0;0;0;0;0;0;5;4;3;2;1;0;5];
nz{7} = [6;5;5;4;4;4;3;3;3;3;2;2;2;2;2;1;1;1;1;1;1;0;0;0;0;0;0;0;6;5;4;3;2;1;0;6];
nz{8} = [7;6;6;5;5;5;4;4;4;4;3;3;3;3;3;2;2;2;2;2;2;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;7;6;5;4;3;2;1;0;7];
nz{9} = [8;7;7;6;6;6;5;5;5;5;4;4;4;4;4;3;3;3;3;3;3;2;2;2;2;2;2;2;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;8;7;6;5;4;3;2;1;0;8];
nz{10} = [9;8;8;7;7;7;6;6;6;6;5;5;5;5;5;4;4;4;4;4;4;3;3;3;3;3;3;3;2;2;2;2;2;2;2;2;1;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;0;9;8;7;6;5;4;3;2;1;0;9];


    for t = 1:NShell_Doublets
        
        mu_begin = Shell_Doublets(t,1);
        mu_end = Shell_Doublets(t,2);
        nu_begin = Shell_Doublets(t,4);
        nu_end = Shell_Doublets(t,5);
        a = Shell_Doublets(t,3);
        b = Shell_Doublets(t,6);
        %If I want to calculate <p|p> or <d|p> type shell doublets, it is
        %much better to use horizontal recursion relations on the already
        %contracted <L1+L2|s> shells. This has to be performed in this
        %program, because the contraction part must be integrated in the
        %HRR part
        basis_a = basis{a};
        basis_b = basis{b};
        
        %{
        if (basis{a}.L > 0 && basis{b}.L > 0)
            for na  = 1:basis{a}.n
                for nb = 1:basis{b}.n
                    %Function recursive_electron_nuclear in electron_nuclear must be overloaded to see
                    %that when g1.L ~ 0 && g2.L ~ 0, it has to generate a
                    %set of Intermediate matrices in order to apply the HRR
                    
                    %------------------------------------------------------
                    %This segment should probably be a function called VRR
                    La = basis{a}.L + basis{b}.L;
                    Lb = 0;
                    
                    g1 = basis{a}.g(na);
                    g2 = basis{b}.g(nb);
                    
                    aa = g1.alpha;
                    ab = g2.alpha;
                    p = aa+ab;
                    P = [aa*g1.x0 + ab*g2.x0,aa*g1.y0 + ab*g2.y0,aa*g1.z0 + ab*g2.z0]/p;
                    C = AL(N,:);
                    RPC2 = sum((C-P).^2);
                    
                    RPA = (P-[g1.x0,g1.y0,g1.z0]);
                    RPB = (P-[g2.x0,g2.y0,g2.z0]);
                    RPC = P-C;
                    
                    Intermediate = zeros((La+1)*(La+2)/2,1);
                    
                    Intermediate = Intermediate + recursive_electron_nuclear(aa,ab,RPA,RPB,RPC,p,RPC2,Z,La,Lb,0);
                    %------------------------------------------------------
                    
                    %Intermediate = Intermediate + electron_nuclear(basis{a}.g(na),basis{b}.g(nb),AL(N,1),AL(N,2),AL(N,3),Z(N))*basis{a}.c(na)*basis{b}.c(nb)*basis{a}.g(na).N*basis{b}.g(nb).N;
                end
            end
            VAB(mu_begin:mu_end,nu_begin:nu_end) = VAB(mu_begin:mu_end,nu_begin:nu_end) + nuc_HRR(Intermediate,basis{a},basis{b});
            
        else
        %}
        if (basis_a.L == 0 && basis_b.L == 0)
            
            for na = 1:basis_a.n
                g1 = basis_a.g(na);
                c1 = basis_a.c(na);
                N1 = g1.N;
                temp = 0;
                for nb = 1:basis_b.n
                    g2 = basis_b.g(nb);
                    c2 = basis_b.c(nb);
                    N2 = g2.N;
                    tempN = 0;
                    for N = 1:nAtom
                        tempN = tempN + electron_nuclear_2(g1,g2,AL(N,1),AL(N,2),AL(N,3),Z(N),Boys_Table);
                    end
                    temp = temp + tempN*c2*N2;
                end
                VAB(mu_begin:mu_end,nu_begin:nu_end) = VAB(mu_begin:mu_end,nu_begin:nu_end) + temp*c1*N1;
            end
        
        elseif (basis_a.L > 0 && basis_b.L > 0)    
            
            %Row vector of distance between nuclei A and B
            Ax = basis_a.g(1).x0;
            Ay = basis_a.g(1).y0;
            Az = basis_a.g(1).z0;
            Bx = basis_b.g(1).x0;
            By = basis_b.g(1).y0;
            Bz = basis_b.g(1).z0;
            RAB = [Bx-Ax,By-Ay,Bz-Az];
            
            if basis_b.L == 1
                switch basis_a.L
                    case 1  %[PP]
                        %TO DO 9 feb 2017
                        %V_ds is contracted from integrals that previously
                        %need uncontracted V_ps, so I'm currently
                        %calculating V_ps two times. It would be better to
                        %only calculate it once.
                        %recursive_electron_nuclear3 already produces two
                        %outputs. It should produce several outputs, and
                        %also use bufferInts.
                        V_ps = Contract_VRR_Nuc(basis_a,basis_b,1,0,AL,Z,Boys_Table); 
                        V_ds = Contract_VRR_Nuc(basis_a,basis_b,2,0,AL,Z,Boys_Table);
                        VAB(mu_begin:mu_end,nu_begin:nu_end) = HRR_Nuc(RAB,V_ds,V_ps,1,1);
                    case 2  %[DP]
                        V_ds = Contract_VRR_Nuc(basis_a,basis_b,2,0,AL,Z,Boys_Table); 
                        V_fs = Contract_VRR_Nuc(basis_a,basis_b,3,0,AL,Z,Boys_Table);
                        VAB(mu_begin:mu_end,nu_begin:nu_end) = HRR_Nuc(RAB,V_fs,V_ds,2,1);
                    case 3  %[FP]
                        V_fs = Contract_VRR_Nuc(basis_a,basis_b,3,0,AL,Z,Boys_Table); 
                  %PROBABLY NOT SUPPORTED!!!
                        V_gs = Contract_VRR_Nuc(basis_a,basis_b,4,0,AL,Z,Boys_Table);
                  %-------------------------      
                        VAB(mu_begin:mu_end,nu_begin:nu_end) = HRR_Nuc(RAB,V_gs,V_fs,3,1);
                end
            elseif basis_b.L == 2
                switch basis_a.L
                    case 1  %[PD]
                    case 2  %[DD]
                    case 3  %[FD]
                end
            elseif basis_b.L == 3
                switch basis_a.L
                    case 1  %[PF]
                    case 2  %[DF]
                    case 3  %[FF]
                end      
            end
                
               
            
 
        else
            L1 = basis_a.L;
            L2 = basis_b.L;
            Dim1 = (L1+1)*(L1+2)/2; %L1 = 1, Dim1 = 3. L1 = 2, Dim1 = 6, L1 = 3, Dim1 = 10, etc
            Dim2 = (L2+1)*(L2+2)/2;
            
            for na  = 1:basis_a.n
                g1 = basis_a.g(na);
                aa = g1.alpha;
                c1 = basis_a.c(na);
                N1 = g1.N;
                temp = zeros(Dim1,Dim2);
                for nb = 1:basis_b.n
                    g2 = basis_b.g(nb);
                    ab = g2.alpha;
                    c2 = basis_b.c(nb);
                    N2 = g2.N;
                    p = aa+ab;
                    Px = (aa*g1.x0 + ab*g2.x0)/p;
                    Py = (aa*g1.y0 + ab*g2.y0)/p;
                    Pz = (aa*g1.z0 + ab*g2.z0)/p;
                    RAB = [g1.x0-g2.x0,g1.y0-g2.y0,g1.z0-g2.z0];
                    RPA = [Px-g1.x0,Py-g1.y0,Pz-g1.z0];
                    RPB = [Px-g2.x0,Py-g2.y0,Pz-g2.z0];
                    tempN = zeros(Dim1,Dim2);
                    for N = 1:nAtom    
                        
                        %P = [aa*g1.x0 + ab*g2.x0,aa*g1.y0 + ab*g2.y0,aa*g1.z0 + ab*g2.z0]/p;
                        %C = [AL(N,1),AL(N,2),AL(N,3)];
                        %RPC2 = sum((C-P).^2);
                        

                        RPC = [Px-AL(N,1),Py-AL(N,2),Pz-AL(N,3)];
                        RPC2 = RPC(1)^2+RPC(2)^2+RPC(3)^2;
                        %Data I need: p, XPA, YPA, ZPA, XPB, YPB, ZPB, RPA2, XPC, YPC, ZPC,Z

                        %[Ex,p,q,Px,Qx] = gprod_1D(g1.x0,g1.alpha,g2.x0,g2.alpha); %Px = x-coordinate Gaussian product center, Qx = x1-x2
                        %[Ey,~,~,Py,Qy] = gprod_1D(g1.y0,g1.alpha,g2.y0,g2.alpha);
                        %[Ez,~,~,Pz,Qz] = gprod_1D(g1.z0,g1.alpha,g2.z0,g2.alpha);

                        %Row vector
%                         RAB = [g1.x0-g2.x0,g1.y0-g2.y0,g1.z0-g2.z0];
%                         RPA = [Px-g1.x0,Py-g1.y0,Pz-g1.z0];
%                         RPB = [Px-g2.x0,Py-g2.y0,Pz-g2.z0];
                        
                        %Order of the integrals. it is zero for the final integral, and
                        %increases for intermediate integrals
                    order = 0;
                    bufferInts = cell(100,1); %The intermediate integrals should be stored here for recursive_electron_nuclear_3
                    Lmax = g1.L+g2.L;
                    tempN = tempN + recursive_electron_nuclear_4(aa,ab,RAB,RPA,RPB,RPC,p,RPC2,Z(N),L1,L2,Lmax,order,Boys_Table,bufferInts);
                    %VAB(mu_begin:mu_end,nu_begin:nu_end) = VAB(mu_begin:mu_end,nu_begin:nu_end) + recursive_electron_nuclear_4(aa,ab,RAB,RPA,RPB,RPC,p,RPC2,Z(N),g1.L,g2.L,Lmax,order,Boys_Table,bufferInts)*c1*c2*N1*N2;
                    end
                    temp = temp + tempN*c2*N2;
                end
                VAB(mu_begin:mu_end,nu_begin:nu_end) = VAB(mu_begin:mu_end,nu_begin:nu_end) + temp*c1*N1;
            end
        end
        %end
    end
end


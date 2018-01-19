function VAB = Contract_VRR_Nuc(basis_a,basis_b,L1,L2,AL,Z,Boys_Table)

%This function contracts Electron_Nuclear attraction matrices
%L1 and L2 are not to be taken from basis_a and basis_b, but given because
%the contract function can be applied to intermediate matrices used for
%horizontal recursion relations. (Eventually they could be used in part of
%the vertical recursion relations).

Dim1 = (L1+1)*(L1+2)/2; %L1 = 1, Dim1 = 3. L1 = 2, Dim1 = 6, L1 = 3, Dim1 = 10, etc
Dim2 = (L2+1)*(L2+2)/2;

nAtom = size(Z,2);

VAB = zeros(Dim1,Dim2);
          
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
                        
                        RPC = [Px-AL(N,1),Py-AL(N,2),Pz-AL(N,3)];
                        RPC2 = RPC(1)^2+RPC(2)^2+RPC(3)^2;
                        
                        %Order of the integrals. it is zero for the final integral, and
                        %increases for intermediate integrals
                        order = 0;
                        %bufferInts = cell(100,1); %The intermediate integrals should be stored here for recursive_electron_nuclear_3
                        Lmax = L1+L2;
                        bufferInts = cell(100,1);
                        tempN = tempN + recursive_electron_nuclear_4(aa,ab,RAB,RPA,RPB,RPC,p,RPC2,Z(N),L1,L2,Lmax,order,Boys_Table,bufferInts);
                    end
                        temp = temp + tempN*c2*N2;
                end
                VAB = VAB + temp*c1*N1;
            end


end
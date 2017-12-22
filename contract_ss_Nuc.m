function VAB = contract_ss_Nuc(basis_a,basis_b,Boys_Table,AL,Z)
nAtom = size(Z,2);
VAB = 0; 
for na = 1:basis_a.n
          g1 = basis_a.g(na);
          c1 = basis_a.c(na);
          aa = g1.alpha;
          N1 = g1.N;
          temp = 0;
      for nb = 1:basis_b.n
              g2 = basis_b.g(nb);
              c2 = basis_b.c(nb);
              ab = g2.alpha;
              N2 = g2.N;
              p = aa + ab;
              
              Px = (aa*g1.x0 + ab*g2.x0)/p;
              Py = (aa*g1.y0 + ab*g2.y0)/p;
              Pz = (aa*g1.z0 + ab*g2.z0)/p;
              RAB = [g1.x0-g2.x0,g1.y0-g2.y0,g1.z0-g2.z0];        
              q = aa*ab/p;
              Kab = exp(-q*(RAB(1)^2+RAB(2)^2+RAB(3)^2));
              
              tempN = 0;
          for N = 1:nAtom
              Cx = AL(N,1);
              Cy = AL(N,2);
              Cz = AL(N,3);
              RPC2 = (Px-Cx)^2+(Py-Cy)^2+(Pz-Cz)^2;

              tempN = tempN - (Z(N)*(2*pi/p)*Kab)*Interpolated_Boys_2(p*RPC2,Boys_Table);
          end
      temp = temp + tempN*c2*N2;
      end
VAB = VAB + temp*c1*N1;
end

end
function [Sab,Tab] = OS1e(KabValues,RPAValues,RPBValues,aValues,bValues,pValues,R12sq)

%This function calculates one electron integrals (overlap and kinetic energy)
%R12sq is a function only on the atoms, not the contracted functions
if (L1 == 0 && L2 == 0) %[s|s]
    Sab = KabValues.*(pi./pValues).^1.5;
    abopValues = aVAlues.*bVAlues./pValues;
    Tab = Sab.*(abopValues.*(3-2*abopValues.*R12sq));
elseif (L1 == 1 && L2 == 0) %[p|s]

    Ss_s = KabValues.*(pi./pValues).^1.5;
    abopValues = aVAlues.*bVAlues./pValues;
    Ts_s= Ss_s.*(abopValues.*(3-2*abopValues.*R12sq));

    Spx_s = PAx.*Ss_s;
    Spy_s = PAy.*Ss_s;
    Spz_s = PAz.*Ss_s;

    Sab = zeros(3,1);
    Sab(1) = sum(Spx_s);
    Sab(2) = sum(Spy_s);
    Sab(3) = sum(Spz_s);

elseif (L1 == 2 && L2 == 0) %[p|s]

    Ss_s = KabValues.*(pi./pValues).^1.5;
    abopValues = aVAlues.*bVAlues./pValues;
    Ts_s= Ss_s.*(abopValues.*(3-2*abopValues.*R12sq));

    Spx_s = PAx.*Ss_s;
    Spy_s = PAy.*Ss_s;
    Spz_s = PAz.*Ss_s;

    S;

    Sab = zeros(6,1);
    Sab(1) = sum(Spx_s);
    Sab(2) = sum(Spy_s);
    Sab(3) = sum(Spz_s);

elseif (L1 == 0 && L2 == 1) %[p|s]

    Ss_s = KabValues.*(pi./pValues).^1.5;
    abopValues = aVAlues.*bVAlues./pValues;
    Ts_s= Ss_s.*(abopValues.*(3-2*abopValues.*R12sq));

    Ss_px = PBx.*Ss_s;
    Ss_py = PBy.*Ss_s;
    Ss_pz = PBz.*Ss_s;

    Sab = zeros(1,3);
    Sab(1) = sum(Ss_px);
    Sab(2) = sum(Ss_py);
    Sab(3) = sum(Ss_pz);
end

end

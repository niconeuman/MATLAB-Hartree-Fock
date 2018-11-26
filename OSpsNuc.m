function Vps = OSpsNuc(RPAValues,RPCValues,VssNValues)

V_px_s_0 = RPAValues(:,1).*VssNValues(:,1)-RPCValues(:,1)*VssNValues(:,2);
V_py_s_0 = RPAValues(:,2).*VssNValues(:,1)-RPCValues(:,2)*VssNValues(:,2);
V_pz_s_0 = RPAValues(:,3).*VssNValues(:,1)-RPCValues(:,3)*VssNValues(:,2);

Vps = zeros(3,1);

Vps(1,1) = sum(V_px_s_0);
Vps(2,1) = sum(V_py_s_0);
Vps(3,1) = sum(V_pz_s_0);


end

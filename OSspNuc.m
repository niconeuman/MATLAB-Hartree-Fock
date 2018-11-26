function Vsp = OSspNuc(RPBValues,RPCValues,VssNValues)

V_s_px_0 = RPBValues(:,1).*VssNValues(:,1)-RPCValues(:,1)*VssNValues(:,2);
V_s_py_0 = RPBValues(:,2).*VssNValues(:,1)-RPCValues(:,2)*VssNValues(:,2);
V_s_pz_0 = RPBValues(:,3).*VssNValues(:,1)-RPCValues(:,3)*VssNValues(:,2);

Vps = zeros(1,3);

Vps(1,1) = sum(V_s_px_0);
Vps(1,2) = sum(V_s_py_0);
Vps(1,3) = sum(V_s_pz_0);


end

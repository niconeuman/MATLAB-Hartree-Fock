function ERI = g_ijkl_driver(mu,nu,kappa,lambda) 

%this function should not work with individual basis functions, but with
%contracted basis functions
global Integral_buffer Contracted_Integral_buffer basis

Integral_Data = [g1.x0,g1.y0,g1.z0,g1.alpha,g1.N,g1.type,...
    g2.x0,g2.y0,g2.z0,g2.alpha,g2.N,g2.type,...
    g3.x0,g3.y0,g3.z0,g3.alpha,g3.N,g3.type,...
    g4.x0,g4.y0,g4.z0,g4.alpha,g4.N,g4.type];

Contracted_Integral_Data = [basis{mu}.type,basis{nu}.type,basis{kappa}.type,basis{lambda}.type];

index = ismember(Contracted_Integral_buffer(:,1:4),Contracted_Integral_Data,'rows');

if Integral_buffer(index,5) ~= 0
    ERI = Integral_buffer(index,5);
else
    if(g1.type == 1 && g2.type == 1 && g3.type == 1 && g4.type == 1)
            ERI = g_SSSS(Integral_Data);
    elseif(g1.type > 1 && g2.type == 1 && g3.type == 1 && g4.type == 1)    
            ERI = g_VRR_1(Integral_Data);
    elseif(g1.type > 1 && g2.type > 1 && g3.type == 1 && g4.type == 1)    
            ERI = g_HRR_2(Integral_Data); %This should 
    elseif(g1.type >= 1 && g2.type >= 1 && g3.type > 1 && g4.type == 1)    
            ERI = g_VRR_1(Integral_Data);
    end
end



end
function [Data,g_ijkl]=g_ijkl_OS(x)

%The program uses Obara-Saika recursion relations to calculate
%ERIs. The idea is that when the program generates a set of integrals, it
%should save new integrals on a buffer matrix, so that next


alpha1 = x(1); %Coeficientes de las exponenciales
alpha2 = x(2);
alpha3 = x(3);
alpha4 = x(4);
type_base1 = x(5); %momento angular
type_base2 = x(6);
type_base3 = x(7);
type_base4 = x(8);
Aug_Bohr_Coordinates1 = [x(9) x(10) x(11)];
Aug_Bohr_Coordinates2 = [x(12) x(13) x(14)];
Aug_Bohr_Coordinates3 = [x(15) x(16) x(17)];
Aug_Bohr_Coordinates4 = [x(18) x(19) x(20)];
R1 = Aug_Bohr_Coordinates1;
R2 = Aug_Bohr_Coordinates2;
R3 = Aug_Bohr_Coordinates3;
R4 = Aug_Bohr_Coordinates4;
if length(x) == 24
d1 = x(21); %coeficientes de contracción
d2 = x(22);
d3 = x(23);
d4 = x(24);
else
    d1 = 1;
    d2 = 1;
    d3 = 1;
    d4 = 1;
end
% [ni,nj|exp(-u^2*r12^2|nk,nl]=Ix*Iy*Iz





%Definiciones generales
%Los índices 1 2 y 3 representan x, y, z
P = (alpha1*R1'+alpha2*R2')/(alpha1+alpha2); %P es el vector posición del centro de gravedad de R1 y R2
Q = (alpha3*R3'+alpha4*R4')/(alpha3+alpha4); %Q es el vector posición del centro de gravedad de R3 y R4

%Función gamma incompleta (orden 0);
S1 = alpha1 + alpha2;
S2 = alpha3 + alpha4;
S4 = S1 + S2;
W = (S1*P+S2*Q)/S4;

%t =
%(alpha1+alpha2)*(alpha3+alpha4)/(alpha1+alpha2+alpha3+alpha4)*norm(P-Q)^2;
%%Idem
t = S1*S2/S4*norm(P-Q)^2;
% du = 0.0001; %esto representa una parte considerable del tiempo de cálculo, cuando las funciones son s.
% u = linspace(0,1,1/du);

%Esto está sacado de la función g_ijkl_factored_numeric, y es para hacer
%una grilla no uniforme (logarítmica)
du = zeros(300,1);

for k = 1:300
    du(k)=0.002*exp((k-1)/25); %80 pts y sigma^2 = 20 suma 20.9 bohr, y parecen razonables, 100 y 25 suma 26.26
end
du = du/sum(du);

u = zeros(length(du),1);
u(1) = du(1)/2;

for k = 2:length(du)
    u(k) = u(k-1) + du(k-1)/2 + du(k)/2;
end

F4 = (u.^(2*4).*exp(-t*u.^2))'*du;
F3 = (2*t*F4+exp(-t))/(2*3+1);
F2 = (2*t*F3+exp(-t))/(2*2+1);
F1 = (2*t*F2+exp(-t))/(2*1+1);
F0 = (2*t*F1+exp(-t))/(2*0+1);



R12 = norm(R2-R1);
R34 = norm(R4-R3);

ap12 = (alpha1*alpha2/(alpha1+alpha2));
ap34 = (alpha3*alpha4/(alpha3+alpha4));

G12 = exp(-ap12*R12^2);
G34 = exp(-ap34*R34^2);

%These normalization factor have to take into account 
Nss1 = ((pi/2/alpha1)^1.5*1/(1))^(-0.5);
Nss2 = ((pi/2/alpha2)^1.5*1/(1))^(-0.5);
Nss3 = ((pi/2/alpha3)^1.5*1/(1))^(-0.5);
Nss4 = ((pi/2/alpha4)^1.5*1/(1))^(-0.5);




function d = kronDel(j,k)
if j == k
    d = 1;
else
    d = 0;
end
end

%Nss1*Nss2*Nss3*Nss4*

%Integrales S, correspondientes a un momento angular total M = 0,1,2,3,4,
%asociadas a las funciones gamma
%Reemplazar zeta por S1
%Reemplazar eta por S2
S0000_0 = d1*d2*d3*d4*G12*G34*2*pi^2.5/(S1*S2*S4^0.5)*F0;
S0000_1 = d1*d2*d3*d4*G12*G34*2*pi^2.5/(S1*S2*S4^0.5)*F1;
S0000_2 = d1*d2*d3*d4*G12*G34*2*pi^2.5/(S1*S2*S4^0.5)*F2;
S0000_3 = d1*d2*d3*d4*G12*G34*2*pi^2.5/(S1*S2*S4^0.5)*F3;
S0000_4 = d1*d2*d3*d4*G12*G34*2*pi^2.5/(S1*S2*S4^0.5)*F4;

if (type_base1 == 0 && type_base2 == 0 && type_base3 == 0 && type_base4 == 0)
%<SS|SS>
    g_ijkl = Nss1*Nss2*Nss3*Nss4*S0000_0;

elseif (type_base1 ~= 0 && type_base2 == 0 && type_base3 == 0 && type_base4 == 0)
    
%<PS|SS> 
Li = type_base1; %reasigno los índices para poder permutarlos si necesito
Lj = type_base2;
Lk = type_base3;
Ll = type_base4;
Size_shell = 3;
g_ijkl = zeros(Size_shell,1);
    
    g_ijkl(1) = (P(1)-R1(1))*S0000_0+(W(1)-P(1))*S0000_1;
    g_ijkl(1) = g_ijkl(1)*Nss1p*Nss1*Nss2*Nss3*Nss4;
    g_ijkl(2) = (P(2)-R1(2))*S0000_0+(W(2)-P(2))*S0000_1;
    g_ijkl(2) = g_ijkl(2)*Nss1p*Nss1*Nss2*Nss3*Nss4;
    g_ijkl(3) = (P(3)-R1(3))*S0000_0+(W(3)-P(3))*S0000_1;
    g_ijkl(3) = g_ijkl(3)*Nss1p*Nss1*Nss2*Nss3*Nss4;

elseif (type_base1 == 0 && type_base2 ~= 0 && type_base3 == 0 && type_base4 == 0)
%<SP|SS>
Li = type_base1; 
Lj = type_base2;
Lk = type_base3;
Ll = type_base4;
Size_shell = 3;
g_ijkl = zeros(Size_shell,1);

    g_ijkl(1) = (P(1)-R2(1))*S0000_0+(W(1)-P(1))*S0000_1;
    g_ijkl(1) = g_ijkl(1)*Nss1p*Nss1*Nss2*Nss3*Nss4;
    g_ijkl(2) = (P(2)-R2(2))*S0000_0+(W(2)-P(2))*S0000_1;
    g_ijkl(2) = g_ijkl(2)*Nss1p*Nss1*Nss2*Nss3*Nss4;
    g_ijkl(3) = (P(3)-R2(3))*S0000_0+(W(3)-P(3))*S0000_1;
    g_ijkl(3) = g_ijkl(3)*Nss1p*Nss1*Nss2*Nss3*Nss4;


elseif (type_base1 == 0 && type_base2 == 0 && type_base3 ~= 0 && type_base4 == 0)
%<SS|PS>

% Li = type_base1; %reasigno los índices para poder permutarlos si necesito
% Lj = type_base2;
% Lk = type_base3;
% Ll = type_base4;
Size_shell = 3;
g_ijkl = zeros(Size_shell,1);

    g_ijkl(1) = (Q(1)-R3(1))*S0000_0+(W(1)-Q(1))*S0000_1;
    g_ijkl(1) = g_ijkl(1)*Nss1p*Nss1*Nss2*Nss3*Nss4;
    g_ijkl(2) = (Q(2)-R3(2))*S0000_0+(W(2)-Q(2))*S0000_1;
    g_ijkl(2) = g_ijkl(2)*Nss1p*Nss1*Nss2*Nss3*Nss4;
    g_ijkl(3) = (Q(3)-R3(3))*S0000_0+(W(3)-Q(3))*S0000_1;
    g_ijkl(3) = g_ijkl(3)*Nss1p*Nss1*Nss2*Nss3*Nss4;

%     g_ijkl = (Q(Li)-R3(Li))*S0000_1+S1/S4*S0000_2;    
%     
 elseif (type_base1 == 0 && type_base2 == 0 && type_base3 == 0 && type_base4 ~= 0)
 %<SS|SP>

% Li = type_base4; %reasigno los índices para poder permutarlos si necesito
% Lj = type_base3;
% Lk = type_base1;
% Ll = type_base2;
Size_shell = 3;
g_ijkl = zeros(Size_shell,1);

    g_ijkl(1) = (Q(1)-R4(1))*S0000_0+(W(1)-Q(1))*S0000_1;
    g_ijkl(1) = g_ijkl(1)*Nss1p*Nss1*Nss2*Nss3*Nss4;
    g_ijkl(2) = (Q(2)-R4(2))*S0000_0+(W(2)-Q(2))*S0000_1;
    g_ijkl(2) = g_ijkl(2)*Nss1p*Nss1*Nss2*Nss3*Nss4;
    g_ijkl(3) = (Q(3)-R4(3))*S0000_0+(W(3)-Q(3))*S0000_1;
    g_ijkl(3) = g_ijkl(3)*Nss1p*Nss1*Nss2*Nss3*Nss4;

%     g_ijkl = (Q(Li)-R4(Li))*S0000_1+S1/S4*S0000_2;
% 
 elseif (type_base1 ~= 0 && type_base2 == 0 && type_base3 ~= 0 && type_base4 == 0)
 %<PS|PS>
% %m = 2
% Li = type_base4; %reasigno los índices para poder permutarlos si necesito
% Lj = type_base3;
% Lk = type_base1;
% Ll = type_base2;

Size_shell = 9;
g_ijkl = zeros(Size_shell,1);
S1010_0 = zeros(Size_shell,1);
%el orden va a ser el siguiente: 
%<PxS|PxS>
%<PxS|PyS>
%<PxS|PzS>
%<PyS|PxS>
%<PyS|PyS>
%<PyS|PzS>
%<PzS|PxS>
%<PzS|PyS>
%<PzS|PzS>

%Una vez que entienda bien esto vectorizarlo, porque es absolutamente
%lineal
    S1000_0(1) = (P(1)-R1(1))*S0000_0+(W(1)-P(1))*S0000_1;
    S1000_0(2) = (P(2)-R1(2))*S0000_0+(W(2)-P(2))*S0000_1;
    S1000_0(3) = (P(3)-R1(3))*S0000_0+(W(3)-P(3))*S0000_1;
    
    %Estos que siguen son 3 componentes, px, py, pz
    S1000_1(1) = (P(1)-R1(1))*S0000_1+(W(1)-P(1))*S0000_2;
    S1000_1(2) = (P(2)-R1(2))*S0000_1+(W(2)-P(2))*S0000_2;
    S1000_1(3) = (P(3)-R1(3))*S0000_1+(W(3)-P(3))*S0000_2;
    
    %De estos hay 9 (listados arriba)
    %<PxS|PS>
    S1010_0(1) = (Q(1)-R3(1))*S1000_1(1)+(W(1)-Q(1))*S1000_1(1)+1/2/S1*(S0000_0-S2/S4*S0000_1);
    S1010_0(2) = (Q(2)-R3(2))*S1000_1(1)+(W(2)-Q(2))*S1000_1(1)+1/2/S1*(S0000_0-S2/S4*S0000_1);
    S1010_0(3) = (Q(3)-R3(3))*S1000_1(1)+(W(3)-Q(3))*S1000_1(1)+1/2/S1*(S0000_0-S2/S4*S0000_1);
    %<PyS|PS>
    S1010_0(4) = (Q(1)-R3(1))*S1000_1(2)+(W(1)-Q(1))*S1000_1(2)+1/2/S1*(S0000_0-S2/S4*S0000_1);
    S1010_0(5) = (Q(2)-R3(2))*S1000_1(2)+(W(2)-Q(2))*S1000_1(2)+1/2/S1*(S0000_0-S2/S4*S0000_1);
    S1010_0(6) = (Q(3)-R3(3))*S1000_1(2)+(W(3)-Q(3))*S1000_1(2)+1/2/S1*(S0000_0-S2/S4*S0000_1);
    %<PzS|PS>
    S1010_0(7) = (Q(1)-R3(1))*S1000_1(3)+(W(1)-Q(1))*S1000_1(3)+1/2/S1*(S0000_0-S2/S4*S0000_1);
    S1010_0(8) = (Q(2)-R3(2))*S1000_1(3)+(W(2)-Q(2))*S1000_1(3)+1/2/S1*(S0000_0-S2/S4*S0000_1);
    S1010_0(9) = (Q(3)-R3(3))*S1000_1(3)+(W(3)-Q(3))*S1000_1(3)+1/2/S1*(S0000_0-S2/S4*S0000_1);   
    
g_ijkl = S1010_0*Nss1p*Nss3p*Nss1*Nss2*Nss3*Nss4;

% %Esta parte es más difícil, porque primero tengo que incrementar el primer
% %índice y luego el segundo, aplicando dos veces la fórmula recursiva, que
% %ahora tendría más términos
%     g_ijkl = (P(Li)-R1(Li))*S0000_2+S2/S4*S0000_3;
% 
% 
% 
% 
% 
%     
%     
%     
%     
%     g_ijkl = La*(S_p_s_12(type_base1)*(S00_34*G00p0(type_base3)+S_p_s_34(type_base3)*G0000)+...
%         S00_12*S00_34*Gp0p0(type_base1,type_base3)+S00_12*S_p_s_34(type_base3)*Gp000(type_base1));
%     
% elseif  (type_base1 == 0 && type_base2 ~= 0 && type_base3 ~= 0 && type_base4 == 0)
% %<SP|PS>    
%     g_ijkl = La*(S_s_p_12(type_base2)*(S00_34*G00p0(type_base3)+S_p_s_34(type_base3)*G0000)+...
%         S00_12*S00_34*G0pp0(type_base2,type_base3)+S00_12*S_p_s_34(type_base3)*G0p00(type_base2));
%     
% elseif (type_base1 ~= 0 && type_base2 == 0 && type_base3 == 0 && type_base4 ~= 0)
% %<PS|SP>    
%     g_ijkl = La*(S_p_s_12(type_base1)*(S00_34*G000p(type_base4)+S_s_p_34(type_base4)*G0000)+...
%         S00_12*S00_34*Gp00p(type_base1,type_base4)+S00_12*S_s_p_34(type_base4)*Gp000(type_base1));    
% 
% elseif (type_base1 == 0 && type_base2 ~= 0 && type_base3 == 0 && type_base4 ~= 0)
% %<SP|SP>    
%     g_ijkl = La*(S_s_p_12(type_base2)*(S00_34*G000p(type_base4)+S_s_p_34(type_base4)*G0000)+...
%         S00_12*S00_34*G0p0p(type_base2,type_base4)+S00_12*S_s_p_34(type_base4)*G0p00(type_base2)); 
%     
% elseif (type_base1 ~= 0 && type_base2 ~= 0 && type_base3 == 0 && type_base4 == 0)
%     
% %<PP|SS>    
%     g_ijkl = La*(S_p_p_12(type_base1,type_base2)*S00_34*G0000+S_p_s_12(type_base1)*S00_34*G0p00(type_base2)+...
%         S_s_p_12(type_base2)*S00_34*Gp000(type_base1)+S00_12*S00_34*Gpp00(type_base1,type_base2));
% 
% elseif (type_base1 == 0 && type_base2 == 0 && type_base3 ~= 0 && type_base4 ~= 0)
% %<SS|PP>    
%     g_ijkl = La*(S_p_p_34(type_base3,type_base4)*S00_12*G0000+S_p_s_34(type_base3)*S00_12*G000p(type_base4)+...
%         S_s_p_34(type_base4)*S00_12*G00p0(type_base3)+S00_12*S00_34*G00pp(type_base3,type_base4));    
%     
%     
% elseif (type_base1 ~= 0 && type_base2 ~= 0 && type_base3 ~= 0 && type_base4 == 0)  
% %<PP|PS>      
%     Gppp0 = -S2/S4^2*(S1*S2/S4*(P(type_base1)-Q(type_base1))*(P(type_base2)-Q(type_base2))*(Q(type_base3)-P(type_base3))*F3+...
%                 1/2*(kronDel(type_base1,type_base2)*(P(type_base3)-Q(type_base3))+...
%                     kronDel(type_base1,type_base3)*(P(type_base2)-Q(type_base2))+...
%                     kronDel(type_base2,type_base3)*(P(type_base1)-Q(type_base1)))*F2);
%     %Gp0pp = Gppp0;
%     
%     g_ijkl = La*(S_p_p_12(type_base1,type_base2)*(S_p_s_34(type_base3)*G0000+S00_34*G00p0(type_base3))+...
%                  S_p_s_12(type_base1)*(S_p_s_34(type_base3)*G0p00(type_base2)+S00_34*G0pp0(type_base2,type_base3))+...
%                  S_s_p_12(type_base2)*(S_p_s_34(type_base3)*Gp000(type_base1)+S00_34*Gp0p0(type_base1,type_base3))+...
%                  S00_12*(S_p_s_34(type_base3)*Gpp00(type_base1,type_base2)+S00_34*Gppp0));
% 
% elseif (type_base1 ~= 0 && type_base2 == 0 && type_base3 ~= 0 && type_base4 ~= 0)
% %<PS|PP>       
%     G0ppp = -S1/S4^2*(S1*S2/S4*(P(type_base3)-Q(type_base3))*(P(type_base4)-Q(type_base4))*(Q(type_base1)-P(type_base1))*F3+...
%                 1/2*(kronDel(type_base3,type_base4)*(P(type_base1)-Q(type_base1))+...
%                     kronDel(type_base3,type_base1)*(P(type_base4)-Q(type_base4))+...
%                     kronDel(type_base4,type_base1)*(P(type_base3)-Q(type_base3)))*F2);
%     Gp0pp = G0ppp;
%     
%     g_ijkl = La*(S_p_p_34(type_base3,type_base4)*(S_p_s_12(type_base1)*G0000+S00_12*Gp000(type_base1))+...
%                  S_p_s_34(type_base3)*(S_p_s_12(type_base1)*G000p(type_base4)+S00_12*G0pp0(type_base4,type_base1))+...
%                  S_s_p_34(type_base4)*(S_p_s_12(type_base1)*G00p0(type_base3)+S00_12*Gp0p0(type_base3,type_base1))+...
%                  S00_34*(S_p_s_12(type_base1)*G00pp(type_base3,type_base4)+S00_12*Gp0pp));             
%              
% elseif (type_base1 == 0 && type_base2 ~= 0 && type_base3 ~= 0 && type_base4 ~= 0)  
% %<SP|PP>       
%     G0ppp = -S1/S4^2*(S1*S2/S4*(P(type_base3)-Q(type_base3))*(P(type_base4)-Q(type_base4))*(Q(type_base2)-P(type_base2))*F3+...
%                 1/2*(kronDel(type_base3,type_base4)*(P(type_base2)-Q(type_base2))+...
%                     kronDel(type_base3,type_base1)*(P(type_base4)-Q(type_base4))+...
%                     kronDel(type_base4,type_base1)*(P(type_base3)-Q(type_base3)))*F2);
%     %Gp0pp = G0ppp;
%     
%     g_ijkl = La*(S_p_p_34(type_base3,type_base4)*(S_s_p_12(type_base2)*G0000+S00_12*Gp000(type_base2))+...
%                  S_p_s_34(type_base3)*(S_s_p_12(type_base2)*G000p(type_base4)+S00_12*G0p0p(type_base4,type_base2))+...
%                  S_s_p_34(type_base4)*(S_s_p_12(type_base2)*G00p0(type_base3)+S00_12*G0pp0(type_base3,type_base2))+...
%                  S00_34*(S_s_p_12(type_base2)*G00pp(type_base3,type_base4)+S00_12*G0ppp));           
% 
% elseif (type_base1 ~= 0 && type_base2 ~= 0 && type_base3 == 0 && type_base4 ~= 0)
% %<PP|SP>    
%         Gppp0 = -S2/S4^2*(S1*S2/S4*(P(type_base1)-Q(type_base1))*(P(type_base2)-Q(type_base2))*(Q(type_base4)-P(type_base4))*F3+...
%                 1/2*(kronDel(type_base1,type_base2)*(P(type_base4)-Q(type_base4))+...
%                     kronDel(type_base1,type_base4)*(P(type_base2)-Q(type_base2))+...
%                     kronDel(type_base2,type_base4)*(P(type_base1)-Q(type_base1)))*F2);
%     %Gp0pp = Gppp0;
%     
%     g_ijkl = La*(S_p_p_12(type_base1,type_base2)*(S_s_p_34(type_base4)*G0000+S00_34*G000p(type_base4))+...
%                  S_p_s_12(type_base1)*(S_s_p_34(type_base4)*G0p00(type_base2)+S00_34*G0p0p(type_base2,type_base4))+...
%                  S_s_p_12(type_base2)*(S_s_p_34(type_base4)*Gp000(type_base1)+S00_34*Gp00p(type_base1,type_base4))+...
%                  S00_12*(S_s_p_34(type_base4)*Gpp00(type_base1,type_base2)+S00_34*Gppp0));
%             
% 
%              
% elseif (type_base1 ~= 0 && type_base2 ~= 0 && type_base3 ~= 0 && type_base4 ~= 0) 
% %<PP|PP>     
%     Gpppp = 1/S4^2*(S1^2*S2^2/S4^2*(P(type_base1)-Q(type_base1))*(P(type_base2)-Q(type_base2))*(Q(type_base3)-P(type_base3))*(Q(type_base4)-P(type_base4))*F4-...
%                     S1*S2/2/S4*(kronDel(type_base1,type_base2)*(P(type_base3)-Q(type_base3))*(P(type_base4)-Q(type_base4))+...
%                     kronDel(type_base1,type_base3)*(P(type_base2)-Q(type_base2))*(P(type_base4)-Q(type_base4))+...
%                     kronDel(type_base1,type_base4)*(P(type_base2)-Q(type_base2))*(P(type_base3)-Q(type_base3))+...
%                     kronDel(type_base2,type_base3)*(P(type_base1)-Q(type_base1))*(P(type_base4)-Q(type_base4))+...
%                     kronDel(type_base2,type_base4)*(P(type_base1)-Q(type_base1))*(P(type_base3)-Q(type_base3))+...
%                     kronDel(type_base3,type_base4)*(P(type_base1)-Q(type_base1))*(P(type_base2)-Q(type_base2)))*F3+...
%                     1/4*(kronDel(type_base1,type_base2)*kronDel(type_base3,type_base4)+...
%                     kronDel(type_base1,type_base3)*kronDel(type_base2,type_base4)+...
%                     kronDel(type_base1,type_base4)*kronDel(type_base2,type_base3))*F2);
%     
%     G0ppp = -S1/S4^2*(S1*S2/S4*(P(type_base1)-Q(type_base1))*(Q(type_base2)-P(type_base2))*(Q(type_base3)-P(type_base3))*F3-...
%                 1/2*(kronDel(type_base1,type_base2)*(P(type_base3)-Q(type_base3))+...
%                     kronDel(type_base1,type_base3)*(P(type_base2)-Q(type_base2))+...
%                     kronDel(type_base2,type_base3)*(P(type_base1)-Q(type_base1)))*F2);
%     Gp0pp = G0ppp;   
%                 
%     Gppp0 = -S2/S4^2*(S1*S2/S4*(P(type_base1)-Q(type_base1))*(P(type_base2)-Q(type_base2))*(Q(type_base3)-P(type_base3))*F3+...
%                 1/2*(kronDel(type_base1,type_base2)*(P(type_base3)-Q(type_base3))+...
%                     kronDel(type_base1,type_base3)*(P(type_base2)-Q(type_base2))+...
%                     kronDel(type_base2,type_base3)*(P(type_base1)-Q(type_base1)))*F2);
%     Gpp0p = Gppp0;
%                 
%                 
%     g_ijkl = La*(S_p_p_12(type_base1,type_base2)*(S_p_p_34(type_base3,type_base4)*G0000+...
%         S_p_s_34(type_base3)*G000p(type_base4)+S_s_p_34(type_base4)*G00p0(type_base3)+S00_34*G00pp(type_base3,type_base4))+...
%         S_p_s_12(type_base1)*(S_p_p_34(type_base3,type_base4)*G0p00(type_base2)+S_p_s_34(type_base3)*G0p0p(type_base2,type_base4)+...
%         S00_34*G0ppp)+S_s_p_12(type_base2)*(S_p_p_34(type_base3,type_base4)*Gp000(type_base1)+S_p_s_34(type_base3)*Gp00p(type_base1,type_base4)+...
%         S_s_p_34(type_base4)*Gp0p0(type_base1,type_base3)+S00_34*Gp0pp)+S00_12*(S_p_p_34(type_base3,type_base4)*Gpp00(type_base1,type_base2)+...
%         S_p_s_34(type_base3)*Gpp0p+S_s_p_34(type_base4)*Gppp0+S00_34*Gpppp));         
% % else
% %     g_ijkl = 0;
end
Data = ones(Size_shell,1)*x(1:20);

end
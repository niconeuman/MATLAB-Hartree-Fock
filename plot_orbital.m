function [Forb,Xg,Yg,Zg] = plot_orbital(g1,Cartesian_type,Npoints,isovalue,eg)
x0 = g1.x0;
y0 = g1.y0;
z0 = g1.z0;
Galpha = g1.alpha;
N = g1.N;
%Cartesian_type is:
%0 s
%1 px
%2 py
%3 pz
%4 dxx
%5 dxy
%6 dxz
%7 dyy
%8 dyz
%9 dzz

Cartesian_powers = [0 0 0;
                    1 0 0;
                    0 1 0;
                    0 0 1;
                    2 0 0;
                    1 1 0;
                    1 0 1;
                    0 2 0;
                    0 1 1;
                    0 0 2];

%The grid should not be uniform, but should have higher density near the region in which the orbital is near the isovalue.
%Maybe the following routine should work
r_mean = sqrt(-log(isovalue)*Galpha); %The 2 should be replaced by some inverse function of the isovalue
%I find this analitically by solving
%isovalue = exp(-Galpha*r^2);
%sqrt(-log(isovalue)*Galpha) = r;
Pow = Cartesian_powers(Cartesian_type+1,:);
 Nx = Npoints*(Pow(1)+1);
 Ny = Npoints*(Pow(2)+1);
 Nz = Npoints*(Pow(3)+1);
% 
% x_mean = r_mean*(Pow(1)+1);
% y_mean = r_mean*(Pow(2)+1);
% z_mean = r_mean*(Pow(3)+1);
% 
% x = [linspace(0,x_mean-1/Galpha,floor(Nx/4)) linspace(x_mean-1/Galpha,x_mean+1/Galpha,Nx)];
% x = [-x x]+x0;
% y = [linspace(0,y_mean-1/Galpha,floor(Ny/4)) linspace(y_mean-1/Galpha,y_mean+1/Galpha,Ny)];
% y = [-y y]+y0;
% z = [linspace(0,z_mean-1/Galpha,floor(Nz/4)) linspace(z_mean-1/Galpha,z_mean+1/Galpha,Nz)];
% z = [-z z]+z0;
%meshgrid in x,y,z
x = linspace(-10/Galpha,10/Galpha,Nx)+x0;
y = linspace(-10/Galpha,10/Galpha,Ny)+y0;
z = linspace(-10/Galpha,10/Galpha,Nz)+z0;

[Xg,Yg,Zg] = meshgrid(x,y,z);


Rsq = (Xg-x0).^2+(Yg-y0).^2+(Zg-z0).^2;
Forb = N*((Xg-x0).^Pow(1).*(Yg-y0).^Pow(2).*(Zg-z0).^Pow(3)).*exp(-Galpha*Rsq);

if strcmp(eg,'x2-y2')
Forbx2 = N*((Xg-x0).^2.*(Yg-y0).^0.*(Zg-z0).^0).*exp(-Galpha*Rsq);    
Forby2 = N*((Xg-x0).^0.*(Yg-y0).^2.*(Zg-z0).^0).*exp(-Galpha*Rsq);
Forb = 1/sqrt(2)*(Forbx2-Forby2);
elseif strcmp(eg,'z2')
Forbx2 = N*((Xg-x0).^2.*(Yg-y0).^0.*(Zg-z0).^0).*exp(-Galpha*Rsq);    
Forby2 = N*((Xg-x0).^0.*(Yg-y0).^2.*(Zg-z0).^0).*exp(-Galpha*Rsq);
Forbz2 = N*((Xg-x0).^0.*(Yg-y0).^0.*(Zg-z0).^2).*exp(-Galpha*Rsq);
Forb = 1/sqrt(2)*(Forbx2+Forby2-2*Forbz2);    
end
%Forb = abs(Forb);


start = 0*ones(3,3);
axes_vectors = 2*eye(3);

quiver3(start(1,1),start(1,2),start(1,3),axes_vectors(1,1),axes_vectors(1,2),axes_vectors(1,3),'Color',[1 0.3 0],'LineWidth',2); hold on;
quiver3(start(2,1),start(2,2),start(2,3),axes_vectors(2,1),axes_vectors(2,2),axes_vectors(2,3),'Color',[0.8 0 0],'LineWidth',2); hold on;
quiver3(start(3,1),start(3,2),start(3,3),axes_vectors(3,1),axes_vectors(3,2),axes_vectors(3,3),'Color',[0.5 0 0],'LineWidth',2); hold on;

%scatter3(0,0,0,500,'filled','MarkerFaceColor','b'); hold on;



%surf(Forb,'EdgeColor','none')
isosurface(Xg,Yg,Zg,Forb,isovalue); hold on
isosurface(Xg,Yg,Zg,Forb,-isovalue);
colormap([0.6,0.8,0;1,0,0])
light
lighting phong
%axis tight equal off
%view(40,30)
camzoom(1.5)
alpha(0.3);
axis([-10 10 -10 10 -10 10]*2/Galpha)

end


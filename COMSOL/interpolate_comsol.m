%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Validation of interpolation method
%%% Import surface=(x,y,z,uy) from Comsol
%%% Interpolate on a regular grid
%%%
%%% Sylvie Jarjayes, 18/09/2018
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

close all;

format long;
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%% Import *.txt file from Comsol: x,y=cste,z,Uy
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

node=round(load('brut.txt'),10); %%%% 1/4 model: 55mm*35mmm, units in m

%%%% Changes coordinates

X=1e3*node(:,1); %%%%%% X=x in mm

Y=1e3*node(:,3); %%%%%%% Y=z in mm

Z=1e6*node(:,4);     %%%%%%% Z=uy vertical displacement in microns
  
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Interpolate: 3 methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

%%%%% Grid for interpolation
 
xlin=linspace(0,55,300);

ylin=linspace(0,35,100);

[x,y]=meshgrid(xlin,ylin);

res=meshgrid(xlin,ylin);

f=scatteredInterpolant(X,Y,Z);


%%% Method to be used


%f.Method='nearest';

%f.Method='linear';

f.Method='natural';

z=f(x,y); %%%%% interpolated surface

%%%% Plot interpolated surface

figure(1);

C=contourf(x,y,z); %%%% contour interpolated surface

colorbar;

figure(2)

surf(x,y,z); %%%% plot interpolated surface in 3D

shading interp;
 
hold on;
  
plot3(X,Y,Z,'mo'); %%%% plot initial points in 3D
% 
% 
% 
% h=gca;
% 
% %h.DataAspectRatio=[20 20 0.1];
% 
colorbar;

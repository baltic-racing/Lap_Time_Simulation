

%Vq = resultierende  eff208 M_eff_inter
%Xq = resultierendes Moment M_M_inter
%Yq = resultierende  Drehzahl M_D_inter

%Define a regular grid and interpolate the scattered data over the grid. 
[M_M_inter,M_D_inter] = meshgrid(0:1:160, 0:1:6000);
M_eff_inter = griddata(M,n,eff208,M_M_inter,M_D_inter);

%Plot the gridded data as a mesh and the scattered data as dots. 
mesh(M_M_inter,M_D_inter,M_eff_inter)
hold on
plot3(M,n,eff208,'o')
xlim([0 160])
ylim([0 6000])


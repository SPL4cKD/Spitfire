clear
clc

W=linspace(0,3e9);             %[Laser intensity (W/m^2)]
R=1; %6671/1.496e8;            %[Radius in AU]
c=3e8;              %[m/s]
dt=10*60;           %[sec]
A=4^2;
m=5:8;               %kg


for i=1:numel(W)
    for j=1:numel(m)
dv(i,j)=(2*W(i))/c*A*dt/m(j)/1000;
    end
end

% rp=1.5e8*1000;
% ra=linspace(rp,1e13);
% a=(rp+ra)/2;
% T=2*pi/sqrt(1.327e20)*a.^1.5/3600/24/365/2;
% 
% plot(ra/1.5e11,T)
% title('Time of Arrival vs distance')
% xlabel('Apoapse Radius (AU)')
% ylabel('Time to Arrival [years]')

plot(W/1e9,dv)
title('\DeltaV vs Laser Power')
xlabel('Laser intensity (GW/m^2)')
ylabel('\DeltaV (km/s)')
legend('5kg','6kg','7kg','8kg')
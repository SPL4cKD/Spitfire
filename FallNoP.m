function [p0, v0, mu, n, scale, cA, bt, graAr, M] = FallNoP

G    =   6.67259e-20; % [km^3/kg-s^2] gravitational constant
M    =   [5.9723, 4.19725e-19, 4.79e-21 ]' * 1e24;
mu   =   M * G;
n    =   size( M, 1 );
scale = [6371.008, 0.1085*1e3, 0.0195*1e4 ];
cA = [4 11 12];

% Time: 2017-SEP-22 00:00:00.000 (Ref Earth)
p0 = [  -6.377850719262221E+03   -5.023336757461595E+01    3.359170072330077E+01; % Earth
        -6.866168454248731E+03   -6.681146443367372E+03    1.374210703989274E+03; % ISS
         2.959732633599232E+04   -4.885719245894279E+04    6.544930814119335E+04];% Chandra

v0 = [   4.335144370399872E-03   -4.267128951345561E-01    1.849765804333802E-01; % Earth
         4.774683631756199E+00    4.272708694755027E-01    6.128856023340462E+00; % ISS
         6.895722891993687E-01    3.619968949574370E-02    2.007321342414172E+00]; % Chandra

pMinus = [-1.080896136799748E+04 -1.695196541324074E+03  3.921740654856121E+02];
vMinus = [8.086137229453791E-03 -4.384005975174894E-01  1.856270724450585E-01];

    
    
p0 = reshape(p0',[1,3*n]);
v0 = reshape(v0',[1,3*n]);

for i = 1:size(p0,2)/3
    p0(1+(i-1)*3:3+(i-1)*3) = p0(1+(i-1)*3:3+(i-1)*3) - pMinus;
    v0(1+(i-1)*3:3+(i-1)*3) = v0(1+(i-1)*3:3+(i-1)*3) - vMinus;
end

p0 = reshape(p0,n*3,1);
v0 = reshape(v0,n*3,1);
graAr = [1, 2, 3];
n = size(graAr,2);
bt = [];
end
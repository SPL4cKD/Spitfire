clc;
clear all;

syms B C E H T W L l z C_1 C_2 H_1 H_2 T_1 T_2 W_1 W_2;

C = C_1 - ((z*C_1)/L) + ((z*C_2)/L);

H = H_1 - ((z*H_1)/L) + ((z*H_2)/L);

T = T_1 - ((z*T_1)/L) + ((z*T_2)/L);

W = W_1 - ((z*W_1)/L) + ((z*W_2)/L);

int(((((E*pi^4)/(32*L^4))*l)*(2*((C*T^3)/12)+(C*T)*((H-T^2)/2)^2+((H/2)-T)^3-(W/3))),z,0,L);
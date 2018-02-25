close all;
clear;
clc;

syms y x P_c P_e

% vpasolve(((2/(y+1))^(1/(y-1))*(G)^(-1/y))/sqrt(((y+1)/(y-1))*(1-(G)^(y-1)/y)) == E, G)
vpasolve(((2/(y+1))^(1/(y-1))*(P_c/P_e)^(1/y))/sqrt((1-(P_c/P_e)^((y-1)/y))) == x, P_c/P_e)
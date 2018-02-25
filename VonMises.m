clear all;
clc;
%% Givens - User Input

prompt = {'Enter sigma 1', 'Enter sigma 2','Enter tau (xy)'};
title = 'Givens';
answer = inputdlg(prompt,title);
sig_psx = str2num(answer{1}); %normal stress in x direction
sig_psy = str2num(answer{2}); %normal stress in y direction 
tau_psxy = str2num(answer{3}); %shear stress on the x face going along the y direction
%% Principal Stresses

%Normal Stresses
sig_ps1 = ((sig_psx + sig_psy)/2)+sqrt(((sig_psx - sig_psy)/2)^2+(tau_psxy)^2) %normal stress 1
sig_ps2 = ((sig_psx + sig_psy)/2)-sqrt(((sig_psx - sig_psy)/2)^2+(tau_psxy)^2) %normal stress 2

%Max Shear
tau_psmax2 = (sig_ps1 - sig_ps2)/2 %maximum shear stress
%% Von Mises Yield Criterion

%Von mises
sig_y = sig_ps1^2 - (sig_ps1*sig_ps2) + sig_ps2^2 + (3*tau_psmax2) %calculation for the von misses 
sig_vmy = sqrt(sig_y) %solution for the von misses
close all
clear
clc


format short
%% SVD Solution to Wahba's Problem
for i=[0.05]                           %Percent Errors Simulated
for j=1:100

A=[0.352 ,  0.864 , 0.360 ; ...
  -0.864 ,  0.152 , 0.480 ; ...
   0.360 , -0.480 , 0.800]  ;          %Rotation Matrix

r1=[ 1.0,    0, 0];
r2=[-0.0,  1.0, 0];
r3=[-0.0, -0.0, 1];
r=[r1;r2;r3];                          %Reference Vectors
n=3;

a=[1,1,1];       %Weights

Max_e=i;
e1=err_gen(Max_e);
e2=err_gen(Max_e);
e3=err_gen(Max_e);

b1=A*r1'+e1';
b2=A*r2'+e2';
b3=A*r3'+e3';

B1=a(1)*b1*r1;
B2=a(2)*b2*r2;
B3=a(3)*b3*r3;

B=B1+B2+B3;

%% Calculating SVD Decomposition

[U,S,V]=svd(B);

d=det(U)*det(V);

Aopt=U*diag([1,1,d])*V;

C=0.5*(abs(b1-Aopt*r1').^2+abs(b2-Aopt*r2').^2+abs(b3-Aopt*r3').^2);
L(j)=sum(C)/numel(C);

end
hold on
histogram(L)
title('100 Monte Carlo Simulations - Orthogonal 3 Vectors')
xlabel('Average Loss Function')
ylabel('# of simulations')
legend('5% Error Magnitude','5% Error Magnitude','1% Error Magnitude')
hold off
Average_loss=sum(L)/numel(L)
end 



















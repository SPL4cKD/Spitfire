close all
clear
clc


format short
%% SVD Solution to Wahba's Problem
for i=[0.05]             %Error Magnitude
for j=1:500

A=[0.352 ,  0.864 , 0.360 ; ...
  -0.864 ,  0.152 , 0.480 ; ...
   0.360 , -0.480 , 0.800]  ;      %Rotation Matrix

r1=[ 1.0,  -0.01, 0];              %Reference Vectors
r2=[ 1.0,   0.01, 0];

a=[1,1,1];                  %Weights

Max_e=i;                    %Generate random error vector of (i) magnitude
e1=err_gen(Max_e);
e2=err_gen(Max_e);

b1=A*r1'+e1';
b2=A*r2'+e2';

B1=a(1)*b1*r1;
B2=a(2)*b2*r2;

B=B1+B2;

%% Calculating SVD Decomposition

[U,S,V]=svd(B);

d=det(U)*det(V);

Aopt=U*diag([1,1,d])*V;

C=0.5*(abs(b1-Aopt*r1').^2+abs(b2-Aopt*r2').^2);
L(j)=sum(C)/numel(C);

end
hold on
histogram(L)
title('100 Monte Carlo Simulations - Co-linear 2 Vectors')
xlabel('Average Loss Function')
ylabel('# of simulations')
legend('5% Error Magnitude','5% Error Magnitude','1% Error Magnitude')
hold off
Average_loss=sum(L)/numel(L)
end 



















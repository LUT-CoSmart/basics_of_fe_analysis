clear all; close all; clc;

E=210000;A=10;L=1000;
F=1000;


K1=KlocRod(A,E,L);
Fvec(1)=0;Fvec(2)=F;


det(K1)

inv(K1)

Kglobred=K1(2,2);

det(Kglobred)



u=(Kglobred)^(-1)*Fvec'
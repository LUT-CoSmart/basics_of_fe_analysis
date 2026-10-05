% Example about numerical integration of
% A(x)=int_0^b x^2 dx
% 
% How many integration points nip are needed to evaluate A(x) accurately? 

clear all;
close all;
clc;
format long;

% number of integration points nip
nip=2;          
b=3;

% gauleg2 computes integration points in -1..1 and weights 
[xiv,wxiv]=gauleg2(-1,1,nip);   

f=0;
for ii=1:nip,
    f=f+wxiv(ii)*fun(xiv(ii),b);
end

f




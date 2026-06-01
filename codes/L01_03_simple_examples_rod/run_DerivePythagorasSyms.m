% How to use MATLAB's symbolic toolbox
% Example: scalars

clear all;
close all;
clc;

syms a b

c=(a^2+b^2)^(1/2);

matlabFunction(c,'file','HypotenusePythagorasFormula','vars',{a,b});


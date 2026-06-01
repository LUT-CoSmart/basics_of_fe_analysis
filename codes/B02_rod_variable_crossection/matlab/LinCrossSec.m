function A = LinCrossSec(A1,A2,L,x)
% Lets compute cross section here 


% Linear interpolation
A = A1*(1 - x/L) + A2*x/L; 
function r=Kinematics_v09(DofsAtNode,q,H,W,L,XI)
% Returns position vector r 
% MKM 20082020: This file is based on old kinem_v004
xi=XI(1);
eta=XI(2);
 
Nvecxi=[(eta*xi)/4 - xi/4 - eta/4 + 1/4;
        xi/4 - eta/4 - (eta*xi)/4 + 1/4;
        eta/4 + xi/4 + (eta*xi)/4 + 1/4;
        eta/4 - xi/4 - (eta*xi)/4 + 1/4];

% Shapefunctions in matrix form
  for ii=1:2
      for jj=1:length(Nvecxi)
          jj2=(jj-1)*2+1+(ii-1);
          Nmxi(ii,jj2)=Nvecxi(jj);
      end
    end

r=Nmxi*q;         % Kinematics is not always defined in this way. Code kinematics for those special elements later, look at an original (old,huge and messy) code written by MKM decade ago.




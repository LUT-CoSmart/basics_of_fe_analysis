function [SS,dEEdu]=Stress2412(uu,X,E,nu,xi,eta)
% Written by MKM at the LUT

dNablaudu=dnablaudu_2412(uu,X,xi,eta);
dNablauTdu=dnablauTdu_2412(uu,X,xi,eta);
nablau=nablau_2412(uu,X,xi,eta);
nablauT=nablau';
% The vector of internal forces based on plane stress
II=diag([1,1]);

% Just for verification
%EE=EE_2412(uu,X,xi,eta);       % Green strain
%dEEdu=dEdu_2412(uu,X,xi,eta);

 
dEEduLin=zeros(2,2,8);  % Lin part
dEEduNLin=zeros(2,2,8); % Nonlin part

for kk=1:8, 
    for ii=1:2,
        for jj=1:2,           
            dEEduLin(ii,jj,kk)=dNablaudu(ii,jj,kk)+dNablauTdu(ii,jj,kk);
                for ll=1:2,
                    dEEduNLin(ii,jj,kk)=dEEduNLin(ii,jj,kk)+(nablauT(ll,ii)*dNablaudu(ll,jj,kk)+dNablauTdu(ll,jj,kk)*nablau(ll,ii));
                end
        end
    end
end

% % Linear strain-displacement relation ##################
%EE=1/2*(nablau+nablau');
%dEEdu=1/2*(dEEduLin);          % lin

% % Nonlinear linear strain-displacement relation ##################
EE=1/2*(nablau+nablauT+nablau'*nablau);
dEEdu=1/2*(dEEduLin+dEEduNLin);
% %#######################################################

eps=[EE(1,1), EE(2,2), 2*EE(1,2)].';

% plain stress assumption sigmazz=0. Suitable for a thin structure.
DD=E/(1-nu^2)*[1, nu, 0;
               nu, 1, 0;
               0, 0, (1-nu)/2];
           
sigma=DD*eps;
SS(1,1)=sigma(1);
SS(2,2)=sigma(2);
SS(1,2)=sigma(3);
SS(2,1)=sigma(3);

% 
% % % plane stress
% lambdaps=E/(1-nu^2);
% SS(1,1)=lambdaps*EE(1,1)+lambdaps*nu*EE(2,2);
% SS(2,2)=lambdaps*nu*EE(1,1)+lambdaps*EE(2,2);
% %SS(1,2)=2*G*EE(1,2);
% %SS(2,1)=2*G*EE(2,1);
% SS(1,2)=lambdaps*((1-nu)/2)*2*EE(1,2);      % sama kuin edellä G:n kanssa
% SS(2,1)=lambdaps*((1-nu)/2)*2*EE(2,1);      % sama kuin edellä G:n kanssa
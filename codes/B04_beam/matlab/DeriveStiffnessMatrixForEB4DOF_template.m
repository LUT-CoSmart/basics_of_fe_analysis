% BK10A6400 Basics of FE Analysis (FEMBasics2025autumn)
% Teacher in charge: Marko Matikainen (MKM)
% Template for the Bonus Task 4.
%
% The code derives stiffness matrix for a two node Eluer Bernoulli beam (4 DOF)element 
% Coded by MKM for student's usage in the FEMBasics2025 course


%  Two node EB element ########################################
%   |-> x 
%   |<---- L ---->|  
%  
%  N1             N2 
%   o-------------o
%   ->            ->
%   u1, theta1    u2, theta2
%
%#######################################################################

% Lets use MATLAB's symbolic toolbox
clear all;
clc;
close all;

% Let's initiliaze variables
syms x L v1 theta1 v2 theta2  A E Iz

% A vector of nodal displacements
uu=[v1, theta1, v2, theta2].';

% Number of nodal displacements of the element
DOFs=4;

% Lecture slides...
% AA= ???;           % TIP: Write constants and uncommet line...you should end 4x4 matrix

% polynomials
%p=[1,x x^2, x^3].';               % Note that Matlab syms wants 
                         % notation .' for transpose. Also command
                         % transpose() works
p=transpose([1,x,x^2,x^3]);
N=p.'*AA^-1;

% Writes function for shapefunctions
matlabFunction(N,'file','Shapef_EB4DOF','vars',{x,L});

vh=N*uu;

% Rotation
theta=diff(vh,x);
% Curvature
Epsxx=diff(theta,x);

% Internal strain energy for a beam element
% Wintdx=????             % TIP: Write energy form Wintdx and uncomment
% line

% Integrate over element's length
Wint=int(Wintdx,x,0,L);


% %#######################################################################
% % Approach 1 K=EI int_0^L d^2 N / dx^2 \cdot d^2 N / dx^2 dx
% for kk=1:DOFs
%     DDN(kk)=diff(diff(N(kk),x),x);
% end
% 
% DDNDDN=DDN.'*DDN;       % this gives you 4x4 matrix. Be careful with syms notations for operators
% 
% for jj=1:DOFs
%     for kk=1:DOFs
%         Kloc(jj,kk)=E*Iz*int(DDNDDN(jj,kk),x,0,L);
%     end
% end
% % %#######################################################################
% % % Approach 1.2: or a bit shorter...
% % for ii=1:DOFs
% %     DDN(ii)=diff(diff(N(ii),x),x);
% %     for jj=1:DOFs
% %         DDNDDN=DDN(ii)*DDN(jj); 
% %         Kloc(ii,jj)=E*Iz*int(DDNDDN,x,0,L);
% %     end
% % end
% % %#######################################################################


%#######################################################################
%#######################################################################
% Approach 2: strain energy -> internal forces -> stiffness matrix

% Vector of internal forces
% Fint= d Wint / d u
for kk=1:DOFs
    Fint(kk)=diff(Wint,uu(kk));
end

% Stiffness matrix
% Kloc= d Fint / d u
for ii=1:DOFs
    for kk=1:DOFs
        Kloc(ii,kk)=diff(Fint(ii),uu(kk));
    end
end

%#######################################################################

% Writes function for the elemental internal energy. 
matlabFunction(Fint,'file','FintEB4DOF','vars',{E,Iz,L,uu});

% Writes function for the elemental stiffness matrix
matlabFunction(Kloc,'file','KlocEB4DOF','vars',{E,Iz,L});

% BK10A6400 Basics of FE Analysis (FEMBasics2022)
% Teacher in charge: Marko Matikainen (MKM)
% Shape functions visualization

% Gives the initial nodal displacements {v=0.1} at nodes 2 and 3
% and approximate displacement {x=0,y=0} at middle of element

clear all
clc;
close all;

% Element size Lx Ly
Lx=2;
Ly=2;

% Parameters for shape functions shown in Cook's book 
a=Lx/2;
b=Ly/2;
%#############################

xvec = -Lx/2:Lx/10:Lx/2;
yvec = -Ly/2:Ly/10:Ly/2;;
xlength=length(xvec);
ylength=length(yvec);

for ii=1:xlength
    for jj=1:ylength
        x=xvec(ii);
        y=yvec(jj);
        xdata(ii,jj)=x;
        ydata(ii,jj)=y;
        %kk=(ii-1)*xlength+jj
         Nvec = ShapefunctionsQ4(Lx,Ly,x,y);
         N1data(ii,jj)=Nvec(1);
         N2data(ii,jj)=Nvec(2);
         N3data(ii,jj)=Nvec(3);
         N4data(ii,jj)=Nvec(4);
        
        % Shape functions form Cook's book Eq (3.4-3)        
        N1dataCook(ii,jj)=(a-x)*(b-y)/(4*a*b);
        N2dataCook(ii,jj)=(a+x)*(b-y)/(4*a*b);
        N3dataCook(ii,jj)=(a+x)*(b+y)/(4*a*b);
        N4dataCook(ii,jj)=(a-x)*(b+y)/(4*a*b);
    end
end


% Let's give displacement {v=0.1} at nodes 2 and 3  
% and approximate displacement {x=0,y=0} at middle of element
uuGivenDispl=[0, 0, 0, 0]';
vvGivenDispl=[0, 0, 1, 1]';

NvecGivenDispl=ShapefunctionsQ4(Lx,Ly,0,0);
uhGivenDispl=NvecGivenDispl*uuGivenDispl
vhGivenDispl=NvecGivenDispl*vvGivenDispl


figure(1)
surf(xdata,ydata,N1data)
hold on
surf(xdata,ydata,N1dataCook)
xlabel('x')
ylabel('y')
title('Shapefunction N_1')

figure(2)
surf(xdata,ydata,N2data)
hold on
surf(xdata,ydata,N2dataCook)
xlabel('x')
ylabel('y')
title('Shapefunction N_2')

figure(3)
surf(xdata,ydata,N3data)
hold on
surf(xdata,ydata,N3dataCook)
xlabel('x')
ylabel('y')
title('Shapefunction N_3')

figure(4)
surf(xdata,ydata,N4data)
hold on
surf(xdata,ydata,N4dataCook)
xlabel('x')
ylabel('y')
title('Shapefunction N_4')
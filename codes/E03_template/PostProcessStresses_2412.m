function PostProcessStresses_2412(P0,uu,xloc,nloc,NumberOfElems,ElemNodes,DofsAtNode,L,H,W,E,nu,DrawingLimits)       
% Function draws stresses at the deformed configuration. (Re)Written by MKM 24082020. 
% The most parts arre based on an original (old,huge and messy) code written by MKM decade ago. 

% Points for drawing.
xiData=[-1,1,1,-1];
etaData=[-1,-1,1,1];
 
%Data=zeros(3,NumberOfElem*4);


 for ii = 1:NumberOfElems
        
        % This is same than in case of residual...avoid dublication, do it smarter later
        X=zeros(ElemNodes*DofsAtNode,1);
        U=zeros(ElemNodes*DofsAtNode,1);
     
         % Elemental nodal displacements for the element ii
        U=uu(xloc(ii,:));      
        % Elemental nodal positions for the element ii from P0
        for jj=1:ElemNodes,
            X((jj-1)*DofsAtNode+1:jj*DofsAtNode)=P0(nloc(ii,jj),:); 
        end         
        X=X(:);
        U=U(:);
        q=X+U;         
           
           
           % Xdata=P0(nloc(ii,:),1);
           % Ydata=P0(nloc(ii,:),2);
            for jj=1:4          % four keypoints in binormalized coordinates -1..1, later more 
                    %kk=(ii-1)*4+jj;
                    xi=xiData(jj);
                    eta=etaData(jj);  
                    
                    r0=Kinematics_v09(DofsAtNode,X,H,W,L,[xi,eta]);    % position at the initial configuration in a physical coordinate system
                    r=Kinematics_v09(DofsAtNode,q,H,W,L,[xi,eta]);     % position at the current configuration in a physical coordinate system 
                    
                    % Returns 2.P-K stress. Uses same functions than in a case of 
                    % elemental residual computation
                    [SS,~]=Stress2412(U,X,E,nu,xi,eta);
                    
                    % Displacement gradient nablau
                    nablau=nablau_2412(U,X,xi,eta);
                    
                    % Deformation gradient
                    II=diag([1,1]);
                    F=II+nablau;
                    J=det(F);
                    % Cauchy stress
                    sigma=1/J*F*SS*F';
                    
                    % Invariants
                    trsigma=trace(sigma);
                    I1=trsigma;
                    I2=1/2*(trsigma^2-trace(sigma^2));
                    
                    % VonMises in 2D
                    sigmaxx=sigma(1,1);sigmayy=sigma(2,2);Tauxy=sigma(1,2);
                    SigmaVM=sqrt(sigmaxx^2+sigmayy^2-sigmaxx*sigmayy+3*Tauxy^2); 
                    
%                     if DIM==3;
%                     I3=det(sigma);
%                     end
                    
                    r0data(jj,:)=r0;
                    Data(jj,1)=r(1);
                    Data(jj,2)=r(2);
                    Data(jj,3)=I1;
                    Data(jj,4)=I2;
                    Data(jj,5)=SigmaVM;
            end           
            
             figure(3)
             %h1=fill(r0data(:,1),r0data(:,2),'blue');   % initial configuration
             %h1.FaceAlpha=0.1;
             h=fill(Data(:,1),Data(:,2),Data(:,3));  
             hold on
             axis equal
%             axis(DrawingLimits)
             set(gca,'XColor', 'none','YColor','none')
            
            
             figure(4)
             %h1=fill(r0data(:,1),r0data(:,2),'blue');   % initial configuration
             %h1.FaceAlpha=0.1;
             h=fill(Data(:,1),Data(:,2),Data(:,4));  
             hold on
             axis equal
 %            axis(DrawingLimits)
             set(gca,'XColor', 'none','YColor','none')
             
             figure(5)
             %h1=fill(r0data(:,1),r0data(:,2),'blue');   % initial configuration
             %h1.FaceAlpha=0.1;
             h=fill(Data(:,1),Data(:,2),Data(:,5));  
             hold on
             axis equal
  %          axis(DrawingLimits)
             set(gca,'XColor', 'none','YColor','none')
        %end        
           
 end
 

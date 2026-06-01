function T = TEB4DOFs(alpha)
% Function produce transformation matrix for a two node rod element

T=[-sin(alpha) cos(alpha) 0 0 0 0;
   0 0 1 0 0 0;
   0 0 0 -sin(alpha) cos(alpha) 0;
   0 0 0 0 0 1];




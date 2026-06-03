function T = TRod(alpha)
% Function produce transformation matrix for a two node rod element

T=[cos(alpha) sin(alpha) 0 0;
   0 0 cos(alpha) sin(alpha)];


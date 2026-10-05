function f=fun(xi,b)

x=b/2*(xi+1);   % scale and translate
f0=x^2;         % original function A(x)
f=f0*b/2;       % dx=b/2 dxi ; no function 
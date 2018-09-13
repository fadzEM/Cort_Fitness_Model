% the reproductive rate function f(rho)
function Y = bx(alpha,B,b,x,x0)
  Y =alpha.* B.*(x<x0) + alpha.* B.*exp(-b.*(x-x0)).*(x>=x0);
 
end

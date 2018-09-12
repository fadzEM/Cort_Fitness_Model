function Y = f(alpha,B,b,x,x0)
% reproductive rate b_x
  Y =alpha.* B.*(x<x0) + alpha.* B.*exp(-b.*(x-x0)).*(x>=x0);
 
end
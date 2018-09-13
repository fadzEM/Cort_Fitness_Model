% The cort vs EC function
function Y = cort(beta, K,n,x)
c_0 = 0.05;
x = 10*x;
y = c_0.*(x<=1) + (c_0+(beta .* x.^n)./(K^n+x.^n)).*(1<x) ;
Y = y;
for i = 1:length(Y)
    if (Y(i) > 1)
       Y(i) = 1;
    end
end    

end
function Y = cort(beta, K,n,x)
% cort as a function of environmental challenges.
x = 10*x;
y = 0.05.*(x<=1) + (0.05+(beta .* x.^n)./(K^n+x.^n)).*(1<x) ;
Y = y;
for i = 1:length(Y)
    if (Y(i) > 1)
       Y(i) = 1;
    end
end    

end
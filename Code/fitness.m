function r = fitness(T, c,c_T, m0, alpha, b, B)
%Intrinsic rate of natural increase
r = (log(B*alpha)./T) - (m0 ./ (1- .5*c)) - b.* (c-c_T)./T;
end
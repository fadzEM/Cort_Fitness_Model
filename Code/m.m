function Y = m(m0,x)
%mortality rate
Y = m0;%./(1-.5.*x) this can be used if the mortality is affected by cort levels
end
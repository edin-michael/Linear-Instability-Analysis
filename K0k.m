function [K0k] = K0k(m,k,h1)
K0k=besselk(m,k.*h1);
end

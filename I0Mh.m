function [I0Mh] = I0Mh(m,M,h)
I0Mh=(besseli(m,(M.*h)));
end

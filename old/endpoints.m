function g = endpoints(f)
%function g = endpoints(f)
%
%ENDPOINTS computes end points of a binary image
%taken from page 507, Digital Image Processing using Matlab by R.C.
%Gonzales, R.E. Woods, and S.L. Eddins
%
% g= ENDPOINTS(F)
persistent lut
if isempty(lut)
    lut = makelut(@endpoint_fcn, 3);
end
g = applylut(f,lut);

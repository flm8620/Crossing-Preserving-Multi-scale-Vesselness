function [b] = bspline(x)
%BSPLINE Summary of this function goes here
%   2nd order centered bspline
x=x+1.5;
b=zeros(size(x));
s0 = x<=0;
s1 = x<=1;
s2 = x<=2;
s3 = x<=3;
c1 = (1-s0).*s1;
c2 = (1-s1).*s2;
c3 = (1-s2).*s3;
b=b+(x.*x).*c1+(-2*(x.*x)+6*x-3).*c2 + (x.*x - 6*x+9).*c3;
b=b*0.5;
end


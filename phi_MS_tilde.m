function [result,fourier] = phi_MS_tilde(x,y,a0,l,N_rho,s_rho,theta_i,N_theta)
cfft2 = @(x) fftshift(fft2(ifftshift(x)));
cifft2 = @(x) fftshift(ifft2(ifftshift(x)));

[X,Y] = meshgrid(x,-y);
[T,R] = cart2pol(X,Y);
T=mod(T,2*pi);
al = a0*exp(l*s_rho);
sx=300*al;
sy=150*al;
theta = 2*pi/N_theta*theta_i;
Rot = [cos(theta), -sin(theta); sin(theta), cos(theta)];
S = Rot*[1/sx^2,0;0,1/sy^2]*Rot';
Gaussian = 1/2/pi/sx/sy*exp(-0.5*(X.^2*S(1,1)+Y.^2*S(2,2)+2*X.*Y*S(1,2)));
Gaussian = Gaussian * al*al;
%figure(2)
%imshow(Gaussian,[])
cX = floor(length(x)/2)+1;
cY = floor(length(y)/2)+1;
A = bspline((mod(T-theta, 2*pi)-pi/2)/s_rho);
A(cX,cY)=1/N_theta;
B0k = bspline(log(R.*al)/s_rho);
if l == N_rho-1
    B0k(cX,cY)=1;
else
    B0k(cX,cY)=0;
end
f = A.*B0k;
phi = cifft2(f);
result = phi/al;
%result = Gaussian.*phi/al;
fourier = cfft2(result);
end


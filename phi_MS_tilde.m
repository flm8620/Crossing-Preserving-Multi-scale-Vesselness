function [result,fourier] = phi_MS_tilde(N,a0,l,N_rho,s_rho,theta_i,N_theta)
cfft2 = @(x) fftshift(fft2(ifftshift(x)));
cifft2 = @(x) fftshift(ifft2(ifftshift(x)));

x = -N:N;
[X,Y] = meshgrid(x,-x);
[T,R] = cart2pol(X,Y);
T=mod(T,2*pi);
al = a0*exp(l*s_rho);
sx=2*N*al;
sy=N*al;
theta = 2*pi/N_theta*theta_i;
Rot = [cos(theta), -sin(theta); sin(theta), cos(theta)];
S = Rot*[1/sx^2,0;0,1/sy^2]*Rot';
Gaussian = 1/2/pi/sx/sy*exp(-0.5*(X.^2*S(1,1)+Y.^2*S(2,2)+2*X.*Y*S(1,2)));
Gaussian = Gaussian * al*al;
%figure(2)
%imshow(Gaussian,[])
cX = floor(length(x)/2)+1;
cY = floor(length(x)/2)+1;
A = bspline((mod(T-theta, 2*pi)-pi/2)/s_rho);
A(cX,cY)=1/N_theta;
B0k = bspline(log(R.*al)/s_rho);
B0k(cX,cY)=0;
%if l == N_rho-1
%    B0k(cX,cY)=0;
%else
%    B0k(cX,cY)=0;
%end
f = A.*B0k;
phi = cifft2(f);
%result = phi/al;
result = Gaussian.*phi/al;
fourier = cfft2(result);
% figure(1)
% imshow(Gaussian,[])
% figure(2)
% imshow(phi,[])
% figure(3)
% imshow(result,[])
% figure(4)
% imshow(fourier,[])
end


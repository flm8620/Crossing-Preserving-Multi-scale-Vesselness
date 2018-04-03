%%
clear all;
cfft2 = @(x) fftshift(fft2(ifftshift(x)));
cifft2 = @(x) fftshift(ifft2(ifftshift(x)));
%%
N_theta = 12;
N_rho = 5;
s_theta = 2*pi/N_theta;
rho = 0:1:500;
a0=0.002;
aN=10;
s_rho = (log(aN)-log(a0))/N_rho;
b=cell(N_rho);
sb = zeros(size(rho));
for l=0:N_rho-1
    al = a0 * exp(l*s_rho);
    b{l+1} = bspline(log(rho.*al)/s_rho);
    plot(rho,b{l+1})
    sb=sb+b{l+1};
    hold on;
end
plot(rho,sb,'DisplayName','sum')
hold off;
%%
N = 512;
x = -N:N;
wavelets_unnorm = zeros(length(x),length(x),N_rho,N_theta);
wavelets = zeros(length(x),length(x),N_rho,N_theta);
F_wavelets_unnorm = zeros(length(x),length(x),N_rho,N_theta);
for l=0:N_rho-1
    for j = 0:N_theta-1
        [phi_MS_t_unnorm, F_phi_MS_t_unnorm] = ...
            phi_MS_tilde(N,a0,l,N_rho,s_rho,j,N_theta);
        wavelets_unnorm(:,:,l+1,j+1)=phi_MS_t_unnorm;
        F_wavelets_unnorm(:,:,l+1,j+1)=F_phi_MS_t_unnorm;
        %imshow(real(F_phi_MS_t_unnorm),[]);
    end
end
F_wavelets_unnorm_sum=squeeze(sum(sum(abs(F_wavelets_unnorm),3),4));
%%
figure(3);
imshow(F_wavelets_unnorm_sum,[]);
% normalization
F_wavelets=F_wavelets_unnorm./F_wavelets_unnorm_sum;
%clearvars F_wavelets_unnorm wavelets_unnorm
for l=0:N_rho-1
    for j = 0:N_theta-1
        wavelets(:,:,l+1,j+1)=cifft2(F_wavelets(:,:,l+1,j+1));
    end
end
%%
for l=0:N_rho-1
    for j = 0:N_theta-1
        imshow(real(wavelets(:,:,l+1,j+1)),[]);
        pause(0.2);
    end
end
%%
I = sum(imread('test_crop2.jpg'),3)/3/255;
scores = zeros(length(x),length(x),N_rho,N_theta);
for l=0:N_rho-1
    for j = 0:N_theta-1
        scores(:,:,l+1,j+1)=conv_fft(wavelets(:,:,l+1,j+1),I);
    end
end
%%
for l=0:N_rho-1
    for j = 0:N_theta-1
        imshow(scores(:,:,l+1,j+1),[]);
        pause(0.2);
    end
end
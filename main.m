%%
clear all;
cfft2 = @(x) fftshift(fft2(ifftshift(x)));
cifft2 = @(x) fftshift(ifft2(ifftshift(x)));
%%
N_theta = 12;
N_rho = 5;
s_theta = 2*pi/N_theta;
rho = 0:1:500;
a0=0.01;
aN=1;
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
N = 128;
x = -N:N;
wavelets = zeros(length(x),length(x),N_rho,N_theta);
F_wavelets = zeros(length(x),length(x),N_rho,N_theta);
for l=0:N_rho-1
    for j = 0:N_theta-1
        [phi_MS_t, F_phi_MS_t] = phi_MS_tilde(x,x,a0,l,N_rho,s_rho,j,N_theta);
        wavelets(:,:,l+1,j+1)=phi_MS_t;
        F_wavelets(:,:,l+1,j+1)=F_phi_MS_t;
        %figure(1)
        %imshow(real(phi_MS_t),[])
        %figure(2)
        %imshow(real(cfft2(phi_MS_t)),[])
    end
end
F_wavelets_sum=squeeze(sum(sum(abs(F_wavelets),3),4));
figure(3);
imshow(F_wavelets_sum,[]);
F_wavelets_sum=repmat(F_wavelets_sum,[1,1,N_rho,N_theta]);
F_wavelets=F_wavelets./F_wavelets_sum;
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
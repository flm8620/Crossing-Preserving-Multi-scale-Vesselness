function [c] = conv_fft(a,b)
cfft2 = @(x) fftshift(fft2(ifftshift(x)));
cifft2 = @(x) fftshift(ifft2(ifftshift(x)));
c=cifft2(cfft2(a).*cfft2(b));
end


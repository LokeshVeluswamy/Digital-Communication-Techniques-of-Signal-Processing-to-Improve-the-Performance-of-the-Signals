clc;
clear all;
close all;
k = 10^7; %Number of symbols
no_taps = 3;
SNR = 0:2:18; 
M=2;
 %Constellation 
for m=1:M 
    H(m,1)=cos(2*pi*(m-1)/M); 
    H(m,2)=sin(2*pi*(m-1)/M); 
end
 %Transmitter
for n = 1:length(SNR)
sym = randi(2,[1,k]);             %Generating N random symbols 
inphase = cos(2*pi*(sym-1)/M);    %In phase Component 
quadrature = sin(2*pi*(sym-1)/M); %Quadrature Component
 
signal = inphase+(1i*(quadrature)); %Composite Signal

f=[0.227, 0.46, 0.688, 0.46, 0.227]; %Channel 2 Taps
sig_with_sym1 = conv(inphase,f);     %Convolution of real part of transmitted symbol & taps of channel
sig_with_sym2 = conv(quadrature,f); %Convolution of real part of transmitted symbol & taps of channel
sig_with_sym = sig_with_sym1+(1i*sig_with_sym2);
% Nosyme with 0db variance 
noise = 1/sqrt(2)*(randn(1,k+length(f)-1) + 1i*randn(1,k+length(f)-1));
 
 % Addition of nosyme to channel
y = sig_with_sym + 10^(-SNR(n)/20)*(noise);
L = length(f);
 
%MMSE equalization 
for taps = 1:no_taps
    tap = 10*taps+1;
hautocorr = conv(f,fliplr(f));
hM = toeplitz([hautocorr([5:end]) zeros(1,2*tap+1-L)], [ hautocorr([5:end]) zeros(1,2*tap+1-L) ]);
hM = hM + 1/2*10^(-SNR(n)/10)*eye(2*tap+1);
d = zeros(1,2*tap+1);
d([-2:2]+tap+1) = fliplr(f);
c_mmse  = [inv(hM)*d.'].';
%Matched filter
yfilt_mmse = conv(y,c_mmse);
yfilt_mmse = yfilt_mmse(tap+2:end);
yfilt_mmse = conv(yfilt_mmse,ones(1,1)); % convolution
filter_sampled_mmse = yfilt_mmse(1:1:k); % sampling at time T
 
 
s_cap=real(filter_sampled_mmse);
s_cap1=imag(filter_sampled_mmse);
 
 
for i = 1:k
if s_cap(i)~= 0 && s_cap(i) > (1/sqrt(2))
s_cap(i) = 1;
else if s_cap(i)~= 0 && s_cap(i) < (-1*(1/sqrt(2)))
s_cap(i) = -1;
else
s_cap(i) = 0;
end
end
 
if s_cap1(i)~= 0 && s_cap1(i) > (1/sqrt(2))
s_cap(i) = 1; 
else if s_cap1(i)~= 0 && s_cap1(i) < (-1*(1/sqrt(2)))
s_cap1(i) = -1;
else
s_cap1(i) = 0;
    end
end
end
s_cap2=s_cap+(1i*(s_cap1)) ;
error(taps,n) = size(find([signal - s_cap2]),2);
 
end
end
 
SER = error/k;                                  %Symbol Error Rate
SER_theoretical = 0.5*erfc(sqrt(10.^(SNR/10))); %Theoretical Error Rate
 
%Plot of SNR and SER
close all;
figure (1)
subplot(2,1,1);
semilogy(SNR,SER_theoretical,'Linewidth',2);
grid on
legend('Theoritical'); 
xlabel('SNR, dB'); 
ylabel('Symbol Error Rate');
title('Symbol error probability for QPSK using MMSE Equalizer Channel 2');
axis tight;
 
subplot(2,1,2);
hold all
semilogy(SNR,SER(1,:),'Linewidth',2);
semilogy(SNR,SER(2,:),'Linewidth',2); 
semilogy(SNR,SER(3,:),'Linewidth',2);
legend('11-Taps','21-Taps','31-Taps');
title('Symbol error probability for BPSK using ZeroForcing Equalizer');
grid on 
axis tight; 

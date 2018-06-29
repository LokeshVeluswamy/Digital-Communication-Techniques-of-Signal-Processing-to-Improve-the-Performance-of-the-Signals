clear all 
clc 
SNR = 0:2:18;
k = 10^7;
M = 2;
f = [0.407, 0.815, 0.407];
m = randi(2,[1,k]);
Es = 1;
inPhase = sqrt(Es)*cos(2*pi*(m-1)/M);
quadrature = sqrt(Es)*sin(2*pi*(m-1)/M);
signal = inPhase+(i*(quadrature));
noisevar = 1/sqrt(2)*(randn(1,k) + i*randn(1,k));
r=conv(signal,f);
% Detection part    
for n = 1:length(SNR)
    rl = r(2:length(r)-1);
    rl = rl + 10^(-SNR(n)/20)*noisevar;
    r1 = real(rl);
    r2 = imag(rl);
    % recieved vectors
    
    %iphat = real(rl)>0;
    
    s_cap = r1;
    s_cap1 = r2;
    for i = 1:k
        if s_cap(i) ~= 0 && s_cap(i) > 1/(sqrt(2))
            s_cap(i) = 1;
        else if s_cap(i)~= 0 && s_cap(i) < (-1/sqrt(2))
            s_cap(i) = -1;
            else
                s_cap(i) = 0;
            end
        end
        
         if s_cap1(i) ~= 0 && s_cap1(i) > 1/(sqrt(2))
            s_cap1(i) = 1;
        else if s_cap1(i)~= 0 && s_cap1(i) < (-1/sqrt(2))
             s_cap1(i) = -1;
            else
                s_cap1(i) = 0;
            end   
         end
    end
    s_cap2 = s_cap + (i*s_cap1);
    error(1,n) = size(find([signal-s_cap2]),2);
    %error1(1,n) = size(find([signal-iphat]),2);
end
SER = error/k;
%SER1 = error1/k;
close all;
hold all
semilogy(SNR,SER,'Linewidth',2);
%semilogy(SNR,SER1,'Linewidth',2);
xlabel('SNR (dB)')
ylabel('symbol error probability')
legend('SER' , 'SER1');
grid on
axis tight;
title('Symbol error probability for given channel');

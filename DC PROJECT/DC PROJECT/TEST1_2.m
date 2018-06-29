clear all 
clc 
SNR = 0:2:18;    %SNR Ratio 
k = 10^7;        % Number of Symbols
M = 2;
f = [0.227, 0.46, 0.688, 0.46, 0.227];      %Channnel 2 Taps
m = randi(2,[1,k]);                         % Random Function Generation
Es = 1;
inPhase = sqrt(Es)*cos(2*pi*(m-1)/M);       % In phase Component Generation
quadrature = sqrt(Es)*sin(2*pi*(m-1)/M);    % Quadrature Component Generation
signal = inPhase+(i*(quadrature));          % Signal Generation
noisevar = 1/sqrt(2)*(randn(1,k) + i*randn(1,k));  % Noise Variance 
r=conv(signal,f);
% Detection Of Signals     
for n = 1:length(SNR)
    rl = r(3:length(r)-2);
    rl = rl + 10^(-SNR(n)/20)*noisevar;
    r1 = real(rl);
    r2 = imag(rl);
    % recieved Signals
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
end
SER = error/k;   %Signal Error Ratio 
close all;       %Plot of SNR to SER
hold all
semilogy(SNR,SER,'Linewidth',2);   %Semilog Plot 
xlabel('SNR (dB)')
ylabel('symbol error probability')
legend('SER')
grid on
axis tight;
title('Symbol error probability for given channel');

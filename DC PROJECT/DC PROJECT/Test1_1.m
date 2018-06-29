clear all 
clc 
SNR = 0:2:18;      %Signal to Noise Ratio 
k = 10^7;          % Number of symbols
M = 2;
f = [0.407, 0.815, 0.407];    %Channel 1 Taps
m = randi(2,[1,k]);           %Random variables
Es = 1;
inPhase = sqrt(Es)*cos(2*pi*(m-1)/M);        %inphase component
quadrature = sqrt(Es)*sin(2*pi*(m-1)/M);     %Quadrature component
signal = inPhase+(i*(quadrature));           %Signal Generation
noisevar = 1/sqrt(2)*(randn(1,k) + i*randn(1,k));   %Noise Variance 
r=conv(signal,f);                             %Convolution of signal 
% Detection of Signals Calculations
for n = 1:length(SNR)
    rl = r(2:length(r)-1);
    rl = rl + 10^(-SNR(n)/20)*noisevar;
    r1 = real(rl);
    r2 = imag(rl);
    % recieved vectors Calculations
    
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
SER = error/k;  %Symbol error ratio
close all;      % Plots of SNR and SER
hold all
semilogy(SNR,SER,'Linewidth',2);    % Semilog Graph
xlabel('SNR (dB)')
ylabel('symbol error probability')
legend('SER' , 'SER1');
grid on
axis tight;
title('Symbol error probability for given channel');

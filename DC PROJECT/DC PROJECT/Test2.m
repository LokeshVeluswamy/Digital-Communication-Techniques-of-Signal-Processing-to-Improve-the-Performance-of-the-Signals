clc;
clear all;
k = 10^7;       %Number of symbols
SNR = 0:2:18;   %Symbol to Noise Ratio
M = 2;
f = [0.407, 0.815, 0.407];   %Channel 1 Taps
no_taps = 3;
 
%  Transmitter Side Calculation
 
for n = 1:length(SNR)
    symbols = randi(2,[1,k]);               %Random symbol generation.
    inphase = cos(2*pi*(symbols-1)/M);      %In Phase Component Generation
    quadrature = sin(2*pi*(symbols-1)/M);   % Quadrature Phase Component
    signal = inphase + 1i*quadrature;       % Signal Generation
    r = conv(signal,f);                     % Convolution Of Signal with Channel Taps
    noise = 1/sqrt(2)*(randn(1,k+length(f)-1) + 1i*randn(1,k+length(f)-1));
    y = r + 10^(-SNR(n)/20)*(noise);        % Adding noise with 0dB
    l = length(f);
    for taps = 1:no_taps
     fm = toeplitz([f([2:end]) zeros(1,10*taps+1-l+1)], [ f([2:-1:1]) zeros(1,10*taps+1-l+1) ]);
     d  = zeros(1,10*taps+1);
     d(taps+1) = 1;
     c = [fm\d.'].';
 
     % Mathched filter
     yfilt = conv(y,c);
     yfilt = yfilt(taps+2:end); 
     yfilt = conv(yfilt,ones(1,1));        % Convolution
     filter_sampled = yfilt(1:1:k);        % Sampling at time T
        
        %Reciever 
        
         s_cap = real(filter_sampled);
         s_cap1 = imag(filter_sampled);
         
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
      %Reciever Hard Decision Coding
      ipHat = real(filter_sampled)>0;
 
      %Counting the Errors
      error1(taps,n) = size(find([signal - ipHat]),2);
      
      
     s_cap2 = s_cap + (i*s_cap1);
     error(taps,n) = size(find([signal-s_cap2]),2);
    end
end
SER = error/k;
SER1 = error1/k;

%Simulated Error Probability
SER_theoretical = 0.5*qfunc(sqrt(10.^(SNR/10)));

%Theoretical_SER
 
% Plot of SNR vs SER
close all 
figure (1)
subplot(3,1,1);
semilogy(SNR,SER_theoretical,'Linewidth',2); 
legend('Theoretical');
title('Symbol error probability for BPSK using ZeroForcing Equalizer');
grid on 
axis tight; 
subplot(3,1,2);
hold all
semilogy(SNR,SER(1,:),'Linewidth',2);
semilogy(SNR,SER(2,:),'Linewidth',2); 
semilogy(SNR,SER(3,:),'Linewidth',2);
legend('11-Taps','21-Taps','31-Taps');
title('Symbol error probability for BPSK using ZeroForcing Equalizer');
grid on 
axis tight;
 
subplot(3,1,3);
hold all
semilogy(SNR,SER1(1,:),'Linewidth',2);
semilogy(SNR,SER1(2,:),'Linewidth',2); 
semilogy(SNR,SER1(3,:),'Linewidth',2); 
 
legend('11-Taps','21-Taps','31-Taps');
title('Symbol error probability for BPSK using ZeroForcing Equalizer hat');
grid on 
axis tight;

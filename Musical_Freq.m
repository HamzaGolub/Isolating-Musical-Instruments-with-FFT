clear all; close all; clc


load 'Clip.mat'

Fs = 44100; %sample rate of the sound clip


S = y'; % transposes in order to have consistent dimensions
w = length(y)/4; % break the spectogram up into four time windows

% isolates the correct windows
S1 = S((1-1)*w+1:1*w); 
S2 = S((2-1)*w+1:2*w); 
S3 = S((3-1)*w+1:3*w); 
S4 = S((4-1)*w+1:4*w); 

L = length(S1)/Fs; % length of each window in seconds
n = length(S1);  % number of elements in each window
t = [0:1/Fs:L - 1/Fs]; % t in sec. relative to the start of the window
tau = 0:0.1:L; % discretization for the Gabor transform
k = 2*pi*(1/L/2)*[0:n/2-1 -n/2:-1]; % discretization in frequency space
ks = fftshift(k); % gotta shift them freqs.

Sgt_spec = zeros(length(ks),length(tau)); % initializing the function for the spectrogram

%Gabor Transform Parameters
a = 400; %this will give you the correct width, so exp(-a(...))
range = [1:1800]; %use this when finding the max and index of the transformed Gabor filtered signal
% i.e., max(TransformedSignal_GaborFiltered(range))

% For S1

for j = 1:length(tau)
   g1 = exp(-a*(t - tau(j)).^2); % Window function
   Sg1 = g1.*S1;
   Sgt1 = fft(Sg1); 
   [M1,ind1] = max(abs(Sgt1(range)));
   filtered = exp(-(1/L)*(abs(k) - k(ind1)).^2).*Sgt1;
   Sgt_spec1(:,j) = fftshift(abs(filtered));
end

% For S2

for j = 1:length(tau)
   g2 = exp(-a*(t - tau(j)).^2); % Window function
   Sg2 = g2.*S2;
   Sgt2 = fft(Sg2); 
   [M2,ind2] = max(abs(Sgt2(range)));
   filtered = exp(-(1/L)*(abs(k) - k(ind2)).^2).*Sgt2;
   Sgt_spec2(:,j) = fftshift(abs(filtered)); 
end

% For S3

for j = 1:length(tau)
   g3 = exp(-a*(t - tau(j)).^2); % Window function
   Sg3 = g3.*S3;
   Sgt3 = fft(Sg3); 
   [M3,ind3] = max(abs(Sgt3(range)));
   filtered = exp(-(1/L)*(abs(k) - k(ind3)).^2).*Sgt3;
   Sgt_spec3(:,j) = fftshift(abs(filtered)); 
end

% For S4

for j = 1:length(tau)
   g4 = exp(-a*(t - tau(j)).^2); % Window function
   Sg4 = g4.*S4;
   Sgt4 = fft(Sg4); 
   [M4,ind4] = max(abs(Sgt4(range)));
   filtered = exp(-(1/L)*(abs(k) - k(ind4)).^2).*Sgt4;
   Sgt_spec4(:,j) = fftshift(abs(filtered)); % We don't want to scale it
end

A1 = Sgt_spec1;   % Shape:  484560x110 double
A2 = Sgt_spec2;    % Shape:  484560x110 double
A3 = Sgt_spec3;    % Shape:  484560x110 double
A4 = Sgt_spec4;    % Shape:  484560x110 double


% Plot of spectrogram for each window

subplot(2,2,1)
pcolor(tau,ks,log10(A1))
shading interp
set(gca,'ylim',[0 800],'Fontsize',16)
colormap(hot)
colorbar
xlabel('time (t)'), ylabel('frequency (k)')
title('S1 Window (Log10 Scale)','Fontsize',16)
 
subplot(2,2,2)
pcolor(tau,ks,log10(A2))
shading interp
set(gca,'ylim',[0 800],'Fontsize',16)
colormap(hot)
colorbar
xlabel('time (t)'), ylabel('frequency (k)')
title('S2 Window (Log10 Scale)','Fontsize',16)
 
subplot(2,2,3)
pcolor(tau,ks,log10(A3))
shading interp
set(gca,'ylim',[0 800],'Fontsize',16)
colormap(hot)
colorbar
xlabel('time (t)'), ylabel('frequency (k)')
title('S3 Window (Log10 Scale)','Fontsize',16)
 
subplot(2,2,4)
pcolor(tau,ks,log10(A4))
shading interp
set(gca,'ylim',[0 800],'Fontsize',16)
colormap(hot)
colorbar
xlabel('time (t)'), ylabel('frequency (k)')
title('S4 Window (Log10 Scale)','Fontsize',16)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Isolate the bassline
S = y'; % transposes in order to have consistent dimensions
L = length(S)/Fs; % total length in sec. 
n = length(S); % total number of elements in S
t = [0:1/Fs:L - 1/Fs]; % time discretization
k = (1/L)*[0:n/2-1 -n/2:-1]; % freq. discretization

% Take the Fourier transform of S, and in freq. space isolate all freqs. 
% (in absolute value) that you determine should be part of the baseline 
% according to spectrogram (or also just by listening); that is, all points
% in the transformed function not within the frequency range you determined
% should be set to zero.

bass = fft(S);

for i = 1:n
    if abs(k(i)) > 200 || abs(k(i)) < 50
        bass(i)=0;
   end
end

% Inverse transform of the thresholded function and saved as A5.

A5 = ifft(bass)';     %Shape:  1938240x1 double


%Plays sound
% p8 = audioplayer(A5, Fs); playblocking(p8);

%Plot of amplitude S over time 

figure(1)
plot(t, A5)
xlabel('time (s)'), ylabel('Amplitude')
title('Amplitude of Isolated Baseline Over Time','Fontsize',16)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Isolate the guitar

S = y'; %reinitialize the S from the previous part above.


guit = fft(S);

for i = 1:n
    if abs(k(i)) > 700 || abs(k(i)) < 300
        guit(i)=0;
   end
end

A6 = ifft(guit)';     %Shape:  1938240x1 double

%Play sound

% p8 = audioplayer(A6, Fs); playblocking(p8);

%Plot of amplitude S over time 

figure(2)
plot(t, A6)
xlabel('time (s)'), ylabel('Amplitude')
title('Amplitude of Isolated Guitar Over Time','Fontsize',16)

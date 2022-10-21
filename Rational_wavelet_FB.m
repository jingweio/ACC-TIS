% Fei Ma
function [t, phi, psi, w, PHI_1, PSI_1, h0, h1, g] = Rational_wavelet_FB(q)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rational Wavelet and FB Construction
% input: q
% output: psi, phi - wavelet and scaling functions; t - time vetor; 
%              w- frequency vector; h0, h1 - lowpass filters, g - highpass        
%              filter; PHI_1, PSI_1 - spectrum of wavelet and scaling functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set fractional scale factor 
p = q+1;
a = p/q;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct frequency specturm 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = -(q+2)*pi:pi/10:(q+2)*pi; 
w_l = length(w);
w_l_half = floor(w_l/2);

% higher space
w_s = w*a; % scale the spectrum by 3/2

for j=1:1:w_l_half
%     PSI_1(j) = sqrt(-1)*Meyer_wavelet(w_s(j), p, q, m_xi);
    PSI_1(j) = Meyer_wavelet(w_s(j), q);
    PHI_1(j) = Meyer_scaling(w_s(j), q);
end

for j = w_l_half+1:w_l
%     PSI_1(j) = -sqrt(-1)*Meyer_wavelet(w_s(j), p, q, m_xi);
    PSI_1(j) = Meyer_wavelet(w_s(j), q);
    PHI_1(j) = Meyer_scaling(w_s(j),  q);
end

% lower space
for j=1:1:w_l_half
%     PSI_2(j) = sqrt(-1)*Meyer_wavelet(w(j), p, q, m_xi);
    PSI_2(j) = Meyer_wavelet(w(j),  q);
    PHI_2(j) = Meyer_scaling(w(j),  q);
end
for j=w_l_half+1:1:w_l
%     PSI_2(j) = -sqrt(-1)*Meyer_wavelet(w(j), p, q, m_xi);
    PSI_2(j) = Meyer_wavelet(w(j),  q);
    PHI_2(j) = Meyer_scaling(w(j),  q);
end

% scaling/wavelet filters construction
H1 = zeros(1, w_l);
H2 = zeros(1, w_l);
G = zeros(1, w_l);

for j = 1:w_l
    if PHI_2(j) ~= 0 
    H1(j) = sqrt(a) * PHI_1(j)/PHI_2(j);
    H2(j) = sqrt(a) * PHI_1(j)/PHI_2(j)*exp(-sqrt(-1)*a*w(j));
    end
    if (PHI_2(j) ~= 0) && (abs(PSI_1(j)) > 1e-10)
    G(j) = sqrt(a) * PSI_1(j)/PHI_2(j);
    end
end

% show frequency spectrum of meyer wavelet and scaling functions
% figure;
% plot(w, abs(PSI_1), '--');
% hold on;
% plot(w, abs(PHI_1));
% title(['magnitude of fractional Meyer wavelet spectrum', '(q = ', int2str(q), ')']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show time domain waveforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = 6;
fs = 3*(q+2);
fs = (q+2)/2;
% fs = (q+2)/3;
 fs = (q+4)/2;
% fs = (q+2)/3*2;
t = -T-0.5:1/(fs*a):T-0.5;
t = -T:1/(fs*a):T;

t_l = length(t);

% sinusoidal bases
y = zeros(1, t_l);
psi = zeros(1, t_l);
phi = zeros(1, t_l);

% sinusoidal bases
y = zeros(1, t_l);
psi = zeros(1, t_l);
phi = zeros(1, t_l);

% psi
for j = 1:1:w_l   
    if w(j) > 0 
        x = PSI_2(j)*exp(i*w(j)*t);
    elseif w(j) < 0
        x = PSI_2(j)*exp(i*w(j)*t);
    else
        x = PSI_2(j)*exp(i*w(j)*t);
    end
    
    y = y + x;
end
psi = real(y);

% phi
y = zeros(1, t_l);
for j = 1:1:w_l
    x = PHI_2(j)*exp(i*w(j)*t);
    y = y + x;
end
phi = real(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% H
y = zeros(1, t_l);
for j = 1:1:w_l
    x = H1(j)*exp(i*w(j)*t);
    y = y + x;
end
h0 = real(y);

y = zeros(1, t_l);
for j = 1:1:w_l
    x = H2(j)*exp(i*w(j)*t);
    y = y + x;
end
h1 = real(y);


% G
y = zeros(1, t_l);
for j = 1:1:w_l
    x = G(j)*exp(i*w(j)*t);
    y = y + x;
end
g = real(y);

% show waveforms
% figure;
% plot(t, psi);
% title(['Meyer\_rational wavelet waveform', '(q = ', int2str(q), ')']);

% figure;
% plot(t, phi);
% title(['Meyer\_rational scaling funcions waveform','(q = ', int2str(q), ')']);


function PSI = Meyer_wavelet(w, q)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate spectrum of meyer wavelet with input frequency value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = q+1;
a = p/q;

w1 = (q-q/(2*q+1))*pi;
w2 = w1*a;
w3 = w2*a;

if (abs(w) >= w1) & (abs(w) <= w2)
    w_mod = q/w1*abs(w) - q;
    V = f_construction(w_mod);
    PSI = (2*pi)^(-1/2) * exp(i*w/2) * sin((pi/2)*V);
else
    if (abs(w) >= w2 ) & (abs(w) <= w3)
        w_mod = q/w2*abs(w) - q;
        V = f_construction(w_mod);
        PSI = (2*pi)^(-1/2) * exp(i*w/2) * cos((pi/2)*V);
    else
        PSI = 0;
    end
end

function PHI = Meyer_scaling(w, q)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate spectrum of meyer wavelet with input frequency value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = q+1;
a = p/q;

w1 = (q-q/(2*q+1))*pi;
w2 = w1*a;
w3 = w2*a;

if (abs(w) <= w1)
    PHI = (2*pi)^(-1/2);
else
    if (abs(w) >= w1 ) & (abs(w) <= w2)
        w_mod = q/w1*abs(w) - q;
        V = f_construction(w_mod);
        PHI = (2*pi)^(-1/2) * cos((pi/2)*V);
    else
        PHI = 0;
    end
end

function V = f_construction(a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construction function (polynomial)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 V = (a^4) * (35 - 84*a + 70*(a^2) - 20*(a^3)); 
% V = 35.*a^4 - 84.*a^5 + 70.*a^6 - 20.*a^7;
% V = 10.*a^3 - 15.*a^4 + 6.*a^5;
% V = 3.*a^2 - 2.*a^3;
% V = (1+ cos(a))/2;

clear;

format long

% -- Laser Properties --
c = 3.*10.^8; % Speed of light in m/s
wavelength = 800.*10.^-9;
k = 2.*pi./wavelength; % Wave vector
freq = 2*pi*(c./wavelength); % Frequency
T = 1./freq; % Period
sigma = 107*10^-15; % FWHM of pulse - altar to match the data and report the final value
%twave1 = .005; % Time offset of wave one to center it

Emax = 100; % Arbitrary

% -- Time Settings --
numsteps = 12*10^4; % Total steps
tf = (10)*sigma; %250*10^-15; % Final time (in seconds)
h = tf./numsteps; % Step size
tvec = zeros(1, numsteps); % Initialize vector
Emax = 1; % Arbitrary
twave1 = 0.5*tf;%initial time of pulse 1 - makes graph look nice
twave2 = 2.*twave1;

for i=2:numsteps % Step through and set time vector
    tvec(i) = tvec(i-1)+h;
end

%-- Tau Settings --
numstepsTau = 12*10^3; %total number of steps
tfTau = tf; %final value of tau
hTau = tfTau./numstepsTau; %step size of tau
tauvec = zeros(1,numstepsTau); %initialize vector

Itot = zeros(1,numstepsTau);

for i=2:hTau % Step through and set tau vector
    tauvec(i) = tauvec(i-1)+htau;
end

i = 8000;

Wave1 = E1(Emax,freq,tvec,twave1,sigma);
Wave2 = E2(Emax,freq,tvec,tauvec(i),twave2,sigma);
E = Wave1 + Wave2;
I = Intensity(E);
line(tvec(1,2:end), E(1,2:end))
%scatter(phase(1,2:end), I(1,2:end))
%scatter(tv(1,2:end), I(1,2:end))
%scatter(tv(1,2:end), E(1,2:end))

% function x = timetopos(t) % Simulates the piezo micrometer, varies 1 micron per 1/2 second
% x = 0.0000000000001*t;
% end

%scatter(tauvec,Itot);

function I = Intensity(E)
I = 4*(E).^2;
end

function E = E2(Emax,freq,tv,tau,twave2,sigma)
steps = length(tv);
E = zeros(1,steps);
for i=1:steps
    E(i) = Emax*cos(freq*(tv(i)-tau-twave2))*exp((-(tv(i)+tau-twave2).^2)./(2*sigma^2));
end
end

function E = E1(Emax,freq,tv,twave1,sigma)
steps = length(tv);
E = zeros(1,steps);
for i=1:steps
    E(i) = Emax*cos(freq*(tv(i)-twave1))*exp(-((tv(i)-twave1).^2)./(2*sigma^2));
end
end
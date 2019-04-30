clear;

format long

% -- Laser Properties --
c = 3.*10.^8; % Speed of light in m/s
wavelength = 800.*10.^-9;
k = 2.*pi./wavelength; % Wave vector
freq = c./wavelength; % Frequency
T = 1./freq; % Period
sigma = .0001; % FWHM of pulse - altar to match the data and report the final value
%twave1 = .005; % Time offset of wave one to center it
twave1 = .0005;
twave2 = .001;

Emax = 10; % Arbitrary

% -- Time Settings --
steps = 30000; % Total steps
%tf = .00000000000005; % Final time
tf = .0015; % Final time
h = tf./steps; % Step size
tv = zeros(1, steps); % Initialize vector

for i=2:steps % Step through and set time vector
    tv(i) = tv(i-1)+h;
end

% -- Tau Settings
tausteps = 300;
tauf = .001;
htau = tauf./tausteps;
tauvec = zeros(1, tausteps);

Itot = zeros(1,tausteps);

for i=2:tausteps % Step through and set tau vector
    tauvec(i) = tauvec(i-1)+htau;
end

i = 100;

Wave1 = E1(Emax,freq,tv,twave1,sigma);
Wave2 = E2(Emax,freq,tv,tauvec(i),twave2,sigma);
E = Wave1 + Wave2;
I = Intensity(E);
scatter(tv(1,2:end), E(1,2:end),1)
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
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
twave1 = 5000.5;

Emax = 1; % Arbitrary

% -- Settings --
steps = 3000000; % Total steps
%tf = .00000000000005; % Final time
tf = .001; % Final time
h = tf./steps; % Step size
tv = zeros(1, steps); % Initialize vector
for i=2:steps % Step through and set ime vector
    tv(i) = tv(i-1)+h;
end
phase = zeros(1, steps); % Initialize phaseshift (position) vector
for i=2:steps % Set the phase vector off of the times using function
    phase(i) = timetopos(tv(i));
end
tau = phase./c;
Wave1 = E1(Emax,freq,tv,twave1);
Wave2 = E2(Emax,freq,tv,tau);
I = Intensity(Wave1,Wave2);

E = Wave1 + Wave2;
%scatter(phase(1,2:end), I(1,2:end))
scatter(tv(1,2:end), I(1,2:end))
%scatter(tv(1,2:end), E(1,2:end))

function x = timetopos(t) % Simulates the piezo micrometer, varies 1 micron per 1/2 second
x = 0.0000000000001*t;
end

function I = Intensity(E1,E2)
I = 4*(E1+E2).^2;
end

function E = E2(Emax,freq,tv,tau)
steps = length(tv);
E = zeros(1,steps);
for i=1:steps
    %E(i) = Emax*cos(k*(0)-freq*tv(i)+ phase(i));
    E(i) = Emax*cos(freq*(tv(i)-tau(i)))*exp((-(tv(i)-tau(i)).^2)./(2*0.0001^2));
end
end

function E = E1(Emax,freq,tv,twave1)
steps = length(tv);
E = zeros(1,steps);
for i=1:steps
    %E(i) = Emax*cos(k*(0)-freq*tv(i));
    %E(i) = Emax*cos(freq*(tv(i)-twave1))*exp(((-tv(i)-twave1).^2)./(2*0.0001^2));
    E(i) = Emax*cos(freq*(tv(i)-twave1))*exp((-(tv(i)-twave1).^2)./(2*0.0001^2));
end
end
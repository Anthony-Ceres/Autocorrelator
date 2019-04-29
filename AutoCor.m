clear;

format long

% -- Laser Properties --
c = 3.*10.^8; % Speed of light in m/s
wavelength = 800.*10.^-9;
k = 2.*pi./wavelength; % Wave vecctor
freq = c./wavelength;
T = 1./freq;

Emax = .1; % Arbitrary

% -- Settings --
Mode = 1; % 1: Both beams, 2: one beam, 3: no beams.
steps = 3000; % Total steps
%tf = .0000000000005; % Final time
tf = 1;
h = tf./steps;
tv = zeros(1, steps);
for i=2:steps
    tv(i) = tv(i-1)+h;
end
phase = zeros(1, steps);
for i=2:steps
    phase(i) = timetopos(tv(i));
end
Wave1 = E1(Emax,k,freq,tv);
Wave2 = E2(Emax,k,freq,tv,phase);
I = Intensity(Wave1,Wave2);

scatter(phase(1,2:end), I(1,2:end))

function x = timetopos(t)
x = 0.0000000000001*sin(2*t);
end

function I = Intensity(E1,E2)
I = 4*(E1+E2).^2;
end

function E = E2(Emax,k,freq,tv,phase)
steps = length(tv);
E = zeros(1,steps);
for i=1:steps
    E(i) = Emax*cos(k*(0)-freq*tv(i)+ phase(i));
end
end

function E = E1(Emax,k,freq,tv)
steps = length(tv);
E = zeros(1,steps);
for i=1:steps
    E(i) = Emax*cos(k*(0)-freq*tv(i));
end
end
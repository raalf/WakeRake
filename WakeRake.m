% Set iterations
iter = 50;                %iterations
count = 0;                %iteration counter
change(1) = 10;           %initial value
cdest = 0.007;            %Cd estimation

% Constants
R = 1717;             % [ft lbf / slug R]
m = 2000;             % Numerical integration divisor

% Airfoil Specific 
tc = 0.1                % thickness ratio of the airfoil 
c = 0.4                 % chord length of airfoil
ksi = 0.15*c              % wakerake distance behind airfoil (recommend 0.15c)

% Pitot tube specific
Do = 0.00262;               % outer pitot tube diameter
Di = 0.0013;              % inner pitot tube diameter
l = .088;                    % length of rake tubes

% Other %
gamma = 1.4;              % ratio of air specific heats
w = 1                 % rake height

% Take Measurements (mean over 5000 samples) %
pt0 =                % freestream total pressure
ptwake =             % avg total pressure in wake
T0 =                 % freestream temperature
q0 =                 % freestream dynamic pressure

% Reduce Data %
p0 = pt0 - q0;              % freestream static pressure
rho0 = p0/(T0*R);         % Freestream density of air
V = sqrt((2*q0)/rho0);    % Freestream velocity
a = sqrt(gamma*R*T0);     % Speed of sound
M = V/a;                  % Mach number
dpt = pt0 - ptwake;         % Average total pressure deficit in wake
cdold = 1;               % Old value of the profile Cd
cdnew = cdest;         % New value of the profile Cd

% Absolute Viscosity Corrected for Temperature
mu = (10^-10*.317)*(T0)^1.5*(734.7/(T0 + 216));

% Chord Reynolds number
Re = (rho0*V*chord)/mu;

% eta calculation, Equation 4.6 from Plaisance
eta_temp = (-0.98e-12*Re^2 + 3.02e-6*Re + 3.74);
eta_temp2 = (sqrt(tc)/sqrt(ksi + 0.3));

% zeta calculation, Equation 4.5 from Plaisance
zeta_temp = (.34e-12*Re^2 - (1.07e-6)*Re + 3.21);
zeta_temp2 = (sqrt(tc)/sqrt(ksi + 0.3));

% (p-p0), Equation 4.9 from Plaisance
dp = (-1.33e-6*Re + 4.36)*((tc)/(0.77 + 3.1*ksi).^2)*q0;

% total pressure deficit corrected for effective displacement of the 
% tube center from equation 2.12 from Plaisance
dptcor_temp = (2*0.131*Do + 2*0.0821*Di)*q0/(w*chord);

i = 1;
while (change(i) >= 0.001 && count < iter)
    count = count + 1;
    cdold = cdnew(i);
    
    % eta calculation, Equation 4.7 from Plaisance
    eta(i) = (-1.08e-12*Re^2 + 3.35e-6*Re + 4.15)*sqrt(cdold*(tc))/(ksi + 0.3);
    
    % zeta calculation, Equation 4.5 from Plaisance
    zeta(i) = zeta_temp.*sqrt(cdold).*zeta_temp2;
    
    % total pressure deficit correction, Equation 2.12 from Plaisance
    dptcor(i) = 2*eta(i)*dptcor_temp;
    dptnew = dpt + dptcor(i);
    
    % trapezoidal rule integration
    dy = zeta(i)/m;
    sum = 0;
    for k = 0:m-1
        y = -0.5*zeta(i) + k*dy;
        inter = (cos(pi*y/(zeta(i))))^2;
        inter2 = (cos(pi*(y + dy)/(zeta(i))))^2;
        sum = sum + 0.5*(sqrt(1 - dp(i)/q0(i) - eta(i)*inter)*(1 - sqrt(1 - eta(i)*inter))+ sqrt(1 - dp(i)/q0(i) - eta(i)*inter2)*(1 - sqrt(1 - eta(i)*inter2)))*dy;
    end
    
    % Integrating factor  eq. 2.11
    F(i) = 4*sum/(eta(i)*zeta(i));
    
    % The new uncorrected value of the profile drag coefficient eq. 2.7
    cdnew(i) = F(i)*w*dptnew/(q0);
    sum2 = 0;
    
    for k = 0:m-1
        y = -0.5*zeta(i) + k*dy;
        inter = (cos(pi*y/(zeta(i))))^2;
        inter2 = (cos(pi*(y+dy)/(zeta(i))))^2;
        sum2 = sum2 + 0.5*(sqrt(1 - dp/q0 - eta(i)*inter)*(1 - sqrt(1 - eta(i)*inter))*...
            (1 + (mach^2)/8*(3*dp/q0 + 3 - 2*gamma - 2*(1 - eta(i)*inter) - (2*gamma - 1)*sqrt(1 - eta(i)*inter))) +...
            sqrt(1 - dp/q0 - eta(i)*inter2)*(1 - sqrt(1 - eta(i)*inter2))*...
            (1 + (mach^2)/8*(3*dp/q0 + 3 - 2*gamma - 2*(1 - eta(i)*inter2) - (2*gamma - 1)*sqrt(1 - eta(i)*inter2))))*dy;
    end
    
    % compressibility correction
    cor = sum2/sum;                            
    cdnew(i) = cdnew(i)*cor;
    
    % change in value through iterations
    change(i) = abs(cdold - cdnew(i))/cdnew(i); 
    
    pt2max(i) = pt0 - eta(i)*q0; %maximum pressure loss in wake
    p2 = p0 + dp; %static pressure in wake
    
    cdimax(i) = sqrt((pt2max(i) - p2)/q0)*(1 - sqrt((pt2max(i) - p0)/q0))/chord;
    
end

c_d(i) = cdnew(i);
if(count == Max_Iter)
    c_d(i) = c_d(i) + 10;
end

% Iterations(i) = count;
% cd*1000-5.7298;
% F-0.8655;
% c_l = ((2 * W) / ( rho0 * V^2 * S));
% c_d;
% F;
% Iterations;
% change;

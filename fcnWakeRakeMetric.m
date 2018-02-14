function [cdnew] = fcnWakeRakeMetric( dpt, T0, q0, rho0 )

%% Tunnel Measurements (From LabJack)
% dpt = .07;             % Average total pressure deficit in wake [Pascal]
% T0 = 293;              % freestream temperature [Kelvin]
% q0 = 5;                % freestream dynamic pressure [Pascal]
% rho0 = 1.225;          % freestream density [kg/m^3]


%% Setting parameters

R = 287.05;           % [J/kgK]
m = 2000;             % Numerical integration divisor



% Airfoil Specific

tc = 0.1;                % thickness ratio of the airfoil
chord = 0.25;            % chord length of airfoil [meter]
ksi = 0.15;              % ratio: distance behind airfoil to chord



% Pitot tube specific
Do = 0.00262;             % outer pitot tube diameter [meter]
Di = 0.0013;              % inner pitot tube diameter [meter]
l = .088;                 % length of rake tubes [meter]



% Other %
gamma = 1.4;          % ratio of air specific heats
w = .15;              % rake height [meter]



%% Reduce Data

V = sqrt((2*q0)/rho0);    % Freestream velocity
a = sqrt(gamma*R*T0);     % Speed of sound
mach = V/a;               % Mach number



% Absolute Viscosity Corrected for Temperature (Sutherland's Law)

mu = 1.716*10^(-5)*(T0/273.15)^(3/2)*((273.15+110.4)/(T0+110.4));



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



%% Preparing for iterations

iter = 50;                %iterations
change(1) = 10;           %initial value
cdest = 0.007;            %Cd estimation



i = 1;

count = 1;                %iteration counter

while (change(i) >= 0.001 && count < iter)
    
    
    
    if count == 1
        
        cdold = cdest;
        
    else
        
        cdold = cdnew(i - 1);
        
    end
    
    
    
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
        
        sum = sum + 0.5*(sqrt(1 - dp/q0 - eta(i)*inter)*...
            (1 - sqrt(1 - eta(i)*inter)) + ...
            sqrt(1 - dp/q0 - eta(i)*inter2)*...
            (1 - sqrt(1 - eta(i)*inter2)))*dy;
        
        
        
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
    change(i + 1) = abs(cdold - cdnew(i))/cdnew(i);
    
    i = i + 1;
    count = count + 1;
    
end

end

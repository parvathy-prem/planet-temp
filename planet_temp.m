% Choose the thermophysical properties you wish to use %

TOP = 0;                % = 1 use "top layer" thermophysical properties from Vasavada et al., 1999
BOTTOM = 0;             % = 1 use "bottom layer" thermophysical properties from Vasavada et al., 1999
TWO_LAYER = 1;          % = 1 Vasavada et al., 1999 two-layer model for thermophysical properties

if TOP==1
    k0 = 9.22e-4;
    x = 1.48;
    rho = 1300;
elseif BOTTOM==1
    k0 = 9.3e-3;
    x = 0.073;
    rho = 1800;
end
c0 = 4200*0.1812;
e = 0.95;
a = 0.11;

% Stefan-Boltzmann constant [W/m2/K4] %
sb = 5.67e-8;

% Solar flux [W/m2] %
Fs = 1366;

% Rotational period [rad/s]; DEFAULT = 2.463e-6;%
w = 0;

% Thermal flux to surface from surroundings [W/m2]; DEFAULT = 250 %
Fsurr = 0;

% Heat flux at lower boundary [W/m2]; DEFAULT = 0.033 %
Qt = 0;

% Spatial resolution [m] %
if TWO_LAYER==0
    % DEFAULT = 0.1*sqrt(k0/(rho*c0*w));
    dz = 0.0001;        
else
    % DEFAULT = 0.1*sqrt(0.5*(9.3e-3+9.22e-4)/(0.5*(1800+1300)*c0*w));
    dz = 0.0001;        
end

% Temporal resolution [s]; DEFAULT = (2*pi/w)/10000;
dt = 0.1;               

% Number of time-steps; DEFAULT = 10000 %
nt = 20000;

% Plot temperature profiles every nt_out timesteps %
nt_out = 2000;

% Number of points in z-direction; DEFAULT = ceil(0.6/dz) %
nz = 200;

% Initial temperature profile %
for i=1:nz
    T_init(i) = 100;
end

T = T_init;
alpha = 0;

% Timestep at which surface enters shadow % 
ntshadow = 100;

% Initialize thermophysical properties %
for i=1:nz
    z(i) = -(i-1)*dz;
    if TWO_LAYER==1
        if z(i)>=-2e-2
            k0 = 9.22e-4;
            x = 1.48;
            rho = 1300;
        else
            k0 = 9.3e-3;
            x = 0.073;
            rho = 1800;
        end
    end
    k(i) = k0*(1 + x*((T(i)/350)^3));
    if T(i)<=350
        f = (T(i)-300)/300;
        c(i) = 4200*(0.1812 +0.1191*f +0.0176*(f^2) +0.2721*(f^3) +0.1869*(f^4));
    else
        c(i) = 4200*(0.2029 +0.0383*(1-exp(-(T(i)-350)/100)));
    end
end

% Initialize arrays to store average, maximum and minimum temperature profiles %
for i=1:nz
    T_avg(i) = 0;
    T_max(i) = 0;
    T_min(i) = 1000;
end

% Temperature profiles are stored in this table %
T_vs_z(:,1) = z.';
col_num = 2;

% Loop over time-steps %
for n=1:nt
    
    % Plot temperature profiles every 2000 time-steps %
    if rem(n,nt_out)==1
        plot(T,z)
        hold on
        T_vs_z(:,col_num) = T.';
        col_num = col_num+1;
    end
    
    % Check if surface has entered shadow %
    if n==ntshadow
        Fs = 0;
    end
    
    % Update time and alpha %
    time(n) = (n-1)*dt;
    Tsurf(n) = T(1);
    
    while alpha>2*pi
        alpha = alpha-2*pi;
    end
    
    alpha = alpha+w*dt;
    
    % START solve heat equation %
    for i=2:nz-1
        if TWO_LAYER==1
            if z(i)>=-2e-2
                k0 = 9.22e-4;
                x = 1.48;
                rho = 1300;
            else
                k0 = 9.3e-3;
                x = 0.073;
                rho = 1800;
            end
        end
        kplus = 0.5*(k(i)+k(i+1));
        kminus = 0.5*(k(i-1)+k(i));
        T_new(i) = T(i) + (dt/(rho*c(i)*dz*dz))*(kplus*(T(i+1)-T(i)) - kminus*(T(i)-T(i-1)));
    end
    
    if TWO_LAYER==1
        k0 = 9.22e-4;
        x = 1.48;
        rho = 1300;
    end
    T_new(1) = T(1) + (0.5*k(2)*dt/(rho*c(1)*dz*dz))*(T(3)-T(1)) - (dt/(rho*c(1)*dz))*(e*sb*(T(1)^4) - ((1-a)*Fs+Fsurr)*max(cos(alpha),0));
    
    if TWO_LAYER==1
        k0 = 9.3e-3;
        x = 0.073;
        rho = 1800;
    end
    T_new(nz) = T(nz) + (0.5*k(nz-1)*dt/(rho*c(nz)*dz*dz))*(T(nz-2)-T(nz)) + (dt/(rho*c(nz)*dz))*Qt;
    
    for i=1:nz
        T(i) = T_new(i);
        if TWO_LAYER==1
            if z(i)>=-2e-2
                k0 = 9.22e-4;
                x = 1.48;
                rho = 1300;
            else
                k0 = 9.3e-3;
                x = 0.073;
                rho = 1800;
            end
        end
        k(i) = k0*(1 + x*((T(i)/350)^3));
        if T(i)<=350
            f = (T(i)-300)/300;
            c(i) = 4200*(0.1812 +0.1191*f +0.0176*(f^2) +0.2721*(f^3) +0.1869*(f^4));
        else
            c(i) = 4200*(0.2029 +0.0383*(1-exp(-(T(i)-350)/100)));
        end
    end
    % END solve heat equation %
    
    % Update average, maximum and minimum temperature profiles %
    for i=1:nz
        T_avg(i) = T_avg(i)+T(i);
        if T(i)>T_max(i)
            T_max(i) = T(i);
        end
        if T(i)<T_min(i)
            T_min(i) = T(i);
        end
    end
    
end

% Calculate average temperature profile %
for i=1:nz
    T_avg(i) = T_avg(i)/nt;
end
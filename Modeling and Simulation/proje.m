clc;close all;clear;

%Initial situation
v0 = 0;                  %initial velocity(m/s)
x0 = 0;                  %initial horizontal position (Range) [meters]
y0 = 0;                  %initial vertical position (Height) [meters]
vx0 = v0;                %horizontal component of velocity(m/s)
vy0 = v0;                %the vertical component of velocity(m/s)
t_burn = 1.5;            %burn time(s) 
m_rocket_initial=10.0;    %rocket mass(kg)
m_propellant = 5;      %propellant mass(kg)
N_monte=2000;            %the number of monte carlo loop  

%Pre-allocation
all_ranges = zeros(N_monte,1);
all_altitudes = zeros(N_monte,1);

%To record a sample trajectory
sample_x = [];
sample_y = [];
sample_t = [];
sample_T = [];
sample_P = [];
sample_rho = [];
sample_h = [];

      
for k=1:1:N_monte

    %Air density
    T_mean = (268.15 + 303.15) / 2;          %mean temperature (Kelvin)
    T_std = 5;                               %standard deviation for temperature (K)
    T0 = T_mean + T_std * randn;             %temperature with normal distribution
    
    P_mean = (895 + 915) / 2 * 100;          %mean pressure (Pa)
    P_std = 2000;                            %standard deviation for pressure (Pa)
    P_0 = P_mean + P_std * randn;            %pressure with normal distribution
    
    R = 287.058;                            %the specific gas constant for dry air
    rho0 = P_0/(R*T0);
    
                               
    
    
    %Drag force constants
    Cd = 0.4;                %aerodynamic drag coefficient   
    r_rocket = 0.078;         %rocket radius(meters)
    A = pi * r_rocket^2;     %Cross-sectional area perpendicular to the direction of flow(m^2) 
     
    %Gravitational acceleration constans           
    G=6.67430 * 10^(-11);           %Universal Gravitational Constant(m^3*kg^-1*s^(-2))
    R_ankara=6.369*10^6;            %Effective Earth radius at Ankara latitude(meters)
    h_alt =938;                     %Average altitude above mean sea level(meters) 
    m_earth=5.9722 * 10^24;         %Earth's Mass(kg)
    r_initial = R_ankara + h_alt;   %Initial radius measured from the center of the Earth
    g0 = (G * m_earth) / (r_initial^2);
    
    %Thrust
    m_dot = m_propellant / t_burn;              %Fuel Flow Rate(kg/s) 
    Isp = 250;                                  %Specific Impulse(s)  
    V_e=Isp*g0;                                 %Exhaust Velocity(m/s)
    P_e = 90000;                                %Nozzle Exit Pressure(Pa)  
    A_e = 0.002;                                %Nozzle Exit Area(m^2)
    
    %Time step
    dt = 0.01;          %time step (s)
    t_max = 150;         %maximum simulation time (s)
    t = 0:dt:t_max;
    N = length(t);
    
    %Allocation
    x = zeros(1,N);
    y = zeros(1,N);
    vx = zeros(1,N);
    vy = zeros(1,N);
    m = zeros(1,N);
    
    % Atmospheric parameters for recording.
    T_array = zeros(1,N);
    P_array = zeros(1,N);
    rho_array = zeros(1,N);
        
    % Initialize the initial conditions
    x(1) = x0;
    y(1) = y0;
    vx(1) = vx0;
    vy(1) = vy0;
    m(1) = m_rocket_initial;
    
    L_lapse = 0.0065;               %Temperature decrease rate (K/m)
    g_expo = g0 / (R * L_lapse);    %Exponent for the pressure formula
    
    
    %%% MAİN SECTİON %%%
    
    for i = 1:N-1
        
        h_current = y(i); 
        
        if y(i) < 0 && i > 1
        break
        end
    
        
        % Update the temperature and pressure change
        if h_current < 11000 
            T_local = T0 - (L_lapse * h_current);
            P_local = P_0 * (1 - (L_lapse * h_current) / T0)^g_expo;
        else
            T_local = 216.65;
            P_local = P_0 * 0.22; 
        end
    
        % Update air density 
        rho_local = P_local / (R * T_local);
        
        % Save
        T_array(i) = T_local;
        P_array(i) = P_local;
        rho_array(i) = rho_local;
   
        % Update mass
        if t(i) <= t_burn
            m(i+1) = m(i) - m_dot * dt;  % Fuel is burning so mass is decreasing
            is_burning = true;
        else
            m(i+1) = m(i);               % Fuel is finished so mass is constant
            is_burning = false;
        end
        
        % Speed vector 
        v_mag = sqrt(vx(i)^2 + vy(i)^2); % Total speed magnitude
    
        %Calculation of forces
    
        %Thrust
        if is_burning
            F_thrust_val = m_dot * V_e + (P_e - P_local) * A_e;  %Total thrust of the engine
        else
            F_thrust_val = 0; 
        end
    
        %Drag force
        F_drag_val = 0.5 * rho_local * Cd * A * v_mag^2;
    
        % Drag force components
        if v_mag ~= 0
            Fdx = -F_drag_val * (vx(i)/v_mag);
            Fdy = -F_drag_val * (vy(i)/v_mag);
        else
            Fdx = 0;
            Fdy = 0;
        end
        
        % Thrust force components
        if is_burning
        theta_launch = 45*pi/180;    
                if v_mag < 0.1
                    Tx = F_thrust_val * cos(theta_launch);
                    Ty = F_thrust_val * sin(theta_launch);
                else
                    Tx = F_thrust_val * (vx(i) / v_mag);
                    Ty = F_thrust_val * (vy(i) / v_mag);
                end
        else
            Tx = 0;
            Ty = 0;
        end
    
        
        % Net forces
        Fx = Tx + Fdx;
        r_current = r_initial + h_current;
        g_current = (G * m_earth) / (r_current^2);
        Fy = Ty + Fdy - m(i)*g_current;
        
    
        
        % Accelerations
        ax = Fx / m(i);
        ay = Fy / m(i);
    
        %Velocity
        vx(i+1) = vx(i) + ax * dt;
        vy(i+1) = vy(i) + ay * dt;
        
        % Location
        x(i+1) = x(i) + vx(i) * dt;
        y(i+1) = y(i) + vy(i) * dt;
    
    end

    % Store results for this Monte Carlo iteration
    all_ranges(k) = max(x);        % Maximum range achieved
    all_altitudes(k) = max(y);     % Maximum altitude achieved

    % Record the data from the first simulation (for the example trajectory)
    if k == 1
       
        valid_indices = 1:i; 
        
        sample_x = x(valid_indices);
        sample_y = y(valid_indices);
        sample_t = t(valid_indices);
        
        sample_T = T_array(valid_indices);
        sample_P = P_array(valid_indices);
        sample_rho = rho_array(valid_indices);
        sample_h = y(valid_indices);
    end
    
end

%% GRAPHS

% Flight Path
figure('Position',[50, 50, 800, 600]);
plot(sample_x, sample_y, 'LineWidth', 2, 'Color', [0.1 0.7 0.2]);
xlabel('Range (m)', 'FontSize', 12);
ylabel('Height (m)', 'FontSize', 12);
title('Rocket Trajectory (Flight Path)', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
axis equal;

% 1. Range for PDF
figure('Position', [100, 100, 800, 600]);
subplot(2,1,1)
histogram(all_ranges, 50, 'Normalization', 'pdf', 'FaceColor', [0.2 0.4 0.8]);
xlabel('Range (m)', 'FontSize', 12);
ylabel('Probability Density', 'FontSize', 12);
title('Range PDF', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% 2. CDF for Range
subplot(2,1,2)
sorted_ranges = sort(all_ranges);
cdf_ranges = (1:N_monte) / N_monte;
plot(sorted_ranges, cdf_ranges, 'LineWidth', 2, 'Color', [0.8 0.2 0.2]);
xlabel('Range (m)', 'FontSize', 12);
ylabel('Cumulative Probability', 'FontSize', 12);
title('Range CDF', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% 3. PDF and CDF for height
figure('Position', [150, 150, 800, 600]);
subplot(2,1,1)
histogram(all_altitudes, 50, 'Normalization', 'pdf', 'FaceColor', [0.2 0.8 0.4]);
xlabel('Height (m)', 'FontSize', 12);
ylabel('Probability Density', 'FontSize', 12);
title('Height PDF', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

subplot(2,1,2)
sorted_altitudes = sort(all_altitudes);
cdf_altitudes = (1:N_monte) / N_monte;
plot(sorted_altitudes, cdf_altitudes, 'LineWidth', 2, 'Color', [0.8 0.4 0.2]);
xlabel('Height (m)', 'FontSize', 12);
ylabel('Cumulative Probability', 'FontSize', 12);
title('Height CDF', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% 4. Time-Height
figure('Position', [200, 200, 800, 600]);
plot(sample_t, sample_y, 'LineWidth', 2, 'Color', [0.1 0.3 0.7]);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Height (m)', 'FontSize', 12);
title('Time - Height Graph', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% 5. Time-Range
figure('Position', [250, 250, 800, 600]);
plot(sample_t, sample_x, 'LineWidth', 2, 'Color', [0.7 0.1 0.3]);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Range (m)', 'FontSize', 12);
title('Time - Range Graph', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% 6. Temperature-Altitude
figure('Position', [300, 300, 800, 600]);
valid_idx = sample_T > 0;
plot(sample_T(valid_idx), sample_h(valid_idx), 'LineWidth', 2, 'Color', [0.9 0.3 0.1]);
xlabel('Temperature (K)', 'FontSize', 12);
ylabel('Height (m)', 'FontSize', 12);
title('Temperature - Height Graph', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% 7. Pressure-Altitude
figure('Position', [350, 350, 800, 600]);
plot(sample_P(valid_idx),sample_h(valid_idx), 'LineWidth', 2, 'Color', [0.3 0.6 0.8]);
xlabel('Temperature (K)', 'FontSize', 12);
xlabel('Pressure (Pa)', 'FontSize', 12);
ylabel('Height (m)', 'FontSize', 12);
title('Pressure - Altitude Graph', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% 8. Temperature-Air Density
figure('Position', [400, 400, 800, 600]);
plot(sample_T(valid_idx), sample_rho(valid_idx), 'LineWidth', 2, 'Color', [0.6 0.2 0.8]);
xlabel('Temperature (K)', 'FontSize', 12);
ylabel('Air Density (kg/m³)', 'FontSize', 12);
title('Temperature - Air Density Graph', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% 9. Statistical Summary
figure('Position', [450, 450, 900, 600]);
subplot(2,2,1)
boxplot(all_ranges);
ylabel('Range (m)', 'FontSize', 11);
title('Range Distribution (Box Plot)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

subplot(2,2,2)
boxplot(all_altitudes);
ylabel('Height (m)', 'FontSize', 11);
title('Height Distribution (Box Plot)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

subplot(2,2,3)
scatter(all_ranges, all_altitudes, 10, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel('Range (m)', 'FontSize', 11);
ylabel('Height (m)', 'FontSize', 11);
title('Range vs. Altitude Relationship', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

subplot(2,2,4)
histogram2(all_ranges, all_altitudes, 30, 'DisplayStyle', 'tile');
xlabel('Range (m)', 'FontSize', 11);
ylabel('Height (m)', 'FontSize', 11);
title('2D Density Map', 'FontSize', 12, 'FontWeight', 'bold');
colorbar;

% Statistical Results
fprintf('\n=== Monte Carlo Simulation Results  ===\n');
fprintf('Number of Simulations: %d\n\n', N_monte);

fprintf('RANGE STATISTICS:\n');
fprintf('Average: %.2f m\n', mean(all_ranges));
fprintf('Standard Deviation: %.2f m\n', std(all_ranges));
fprintf('Minimum: %.2f m\n', min(all_ranges));
fprintf('Maximum: %.2f m\n', max(all_ranges));
fprintf('Median: %.2f m\n', median(all_ranges));

fprintf('\nHEIGHT STATISTICS:\n');
fprintf('Average: %.2f m\n', mean(all_altitudes));
fprintf('Standard Deviation: %.2f m\n', std(all_altitudes));
fprintf('Minimum: %.2f m\n', min(all_altitudes));
fprintf('Maximum: %.2f m\n', max(all_altitudes));
fprintf('Median: %.2f m\n', median(all_altitudes));

%% initialize
N = 10000;
h = 0.1;
U = 5;

% defines for Nomoto
T = 20;
k = 0.1;
bias = 0.001;

delta_numerator = [k*U];
delta_denominator = [T 1 0 0];
h1 = tf(delta_numerator, delta_denominator);

bias_numerator = [U];
bias_denominator = [T 1 0 0];
h2 = tf(bias_numerator, bias_denominator);

% initial conditions
x0 = 0; % meter
y0 = 100; % meter

psi0 = 0; % rad
r0 = 0; % rad/s

% define functions (given assumptions)
x_dot_func = @(psi, U) [U*cos(psi)];
y_dot_func = @(psi, U) [U*sin(psi)];

% states
x = zeros(N+1, 1);
y = zeros(N+1, 1);
x_dot = zeros(N+1, 1);
y_dot = zeros(N+1, 1);
t = zeros(N+1, 1);
psi = zeros(N+1, 1);
delta = zeros(N+1, 1);
r = zeros(N+1, 1);
r_dot = zeros(N+1, 1);

% pid-variables
k_p = 0.0005;
k_i = 0.000001;
k_d = 0.03;

% state init
x(1) = x0;
y(1) = y0;
x_dot(1) = 0;
y_dot(1) = 0;
psi(1) = 0;
r(1) = 0;
r_dot(1) = 0;

for i = 1:N+1
    % time starts at 0, then increment by the timestep
    t(i) = (i-1) * h; 
    
    % pid-controller
    delta_init = -k_p*y(i)-k_d*U*psi(i)-k_i*trapz(y);
    delta(i) = min([abs(delta_init) 20]) * sign(delta_init);
    
    % simulate 
    r_dot(i) = k*delta(i)/T + bias/T - r(i)/T;
    r(i) = euler2(r_dot(i), r(i), h);
    
    psi(i) = euler2(r(i), psi(i), h); 
    
    x_dot(i+1) = x_dot_func(psi(i), U);
    y_dot(i+1) = y_dot_func(psi(i), U);
    
    % integrate to get new states
    x(i+1) = euler2(x_dot(i+1),x(i),h); 
    y(i+1) = euler2(y_dot(i+1),y(i),h); 
end

%% PLOTTING

x = x(1:end-1);
y = y(1:end-1);
x_dot = x_dot(1:end-1);
y_dot = y_dot(1:end-1);

figure (1); clf;
plot(t, psi, 'b');
grid on;
legend('\psi');
title('Yaw angle');
xlabel('time [s]'); 
ylabel('Degrees');

figure (2); clf;
hold on;
plot(t, delta, 'b');
hold off;
grid on;
legend('\delta');
title('Input');
xlabel('time [s]'); 
ylabel('Degrees');

figure (3); clf;
hold on;
plot(t, r, 'b');
hold off;
grid on;
legend('$\dot{\psi}$','Interpreter','latex');
title('Yaw rate');
xlabel('time [s]'); 
ylabel('Degrees/s');

figure (100); clf;
hold on;
plot(x, y, 'b');
hold off;
grid on;
legend('x');
title('Position plot');
xlabel('x [m]'); 
ylabel('y [m]');



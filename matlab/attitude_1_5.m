% M-script for numerical integration of the attitude dynamics of a rigid 
% body represented by unit quaternions. The MSS m-files must be on your
% Matlab path in order to run the script.
%
% System:                      .
%                              q = T(q)w
%                              .
%                            I w - S(Iw)w = tau
% Control law:
%                            tau = constant
% 
% Definitions:             
%                            I = inertia matrix (3x3)
%                            S(w) = skew-symmetric matrix (3x3)
%                            T(q) = transformation matrix (4x3)
%                            tau = control input (3x1)
%                            w = angular velocity vector (3x1)
%                            q = unit quaternion vector (4x1)
%
% Author:                   2018-08-15 Thor I. Fossen and H�kon H. Helgesen

%% USER INPUTS
h = 0.1;                     % sample time (s)
N  = 4500;                   % number of samples. Should be adjusted

% model parameters
m = 180;
r = 2;
I = m*r^2*eye(3);            % inertia matrix
I_inv = inv(I);

% constants
deg2rad = pi/180;   
rad2deg = 180/pi;

phi = -10*deg2rad;            % initial Euler angles
theta = 10*deg2rad;
psi = 5*deg2rad;

q = euler2q(phi,theta,psi);   % transform initial Euler angles to q

w = [0 0 0]';                 % initial angular rates

table = zeros(N+1,17);        % memory allocation

k_d = 400;
k_p = 20;
K_d = k_d * eye(3);

%% SIMULATION LOOP
for i = 1:N+1
   t = (i-1)*h;   % time
   q_d = euler2q(0,15*cos(0.1*t)*deg2rad, 10*sin(0.05*t)*deg2rad); % reference
   
   q_squiggle = quatprod(quatconj(q_d), q);      % reference error
   tau = -K_d*w - k_p*q_squiggle(2:4); % control law from 1.5

   [phi,theta,psi] = q2euler(q); % transform q to Euler angles
   [phi_squiggle,theta_squiggle,psi_squiggle] = q2euler(q_d); % transform q to Euler angles
   [J,J1,J2] = quatern(q);       % kinematic transformation matrices
   
   q_dot = J2*w;                        % quaternion kinematics
   w_dot = I_inv*(Smtrx(I*w)*w + tau);  % rigid-body kinetics
   
   table(i,:) = [t q' phi theta psi w' tau' phi_squiggle theta_squiggle psi_squiggle];  % store data in table
   
   q = q + h*q_dot;	             % Euler integration
   w = w + h*w_dot;
   
   q  = q/norm(q);               % unit quaternion normalization
end

%% PLOT FIGURES
t       = table(:,1);  
q       = table(:,2:5); 
phi     = rad2deg*table(:,6);
theta   = rad2deg*table(:,7);
psi     = rad2deg*table(:,8);
w       = rad2deg*table(:,9:11);  
tau     = table(:,12:14);
phi_squiggle   = rad2deg*table(:,15);
theta_squiggle = rad2deg*table(:,16);
psi_squiggle   = rad2deg*table(:,17);

phi_d   = 0;
theta_d = 15*cos(0.1*t);
psi_d   = 10*sin(0.05*t);


fig1 = figure (1); clf;
hold on;
% states
plot(t, phi, 'b');
plot(t, theta, 'r');
plot(t, psi, 'g');
%references
plot(t, phi_d, '--b');
plot(t, theta_d, '--r');
plot(t, psi_d, '--g');

hold off;
grid on;
legend('\phi', '\theta', '\psi');
title('Euler angles');
xlabel('time [s]'); 
ylabel('angle [deg]');

set(fig1, 'Units', 'Inches');
pos1 = get(fig1, 'Position');
set(fig1, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos1(3), pos1(4)]);
print(fig1, '1_5_euler_angles', '-dpdf', '-r0');


fig2 = figure (2); clf;
hold on;
plot(t, w(:,1), 'b');
plot(t, w(:,2), 'r');
plot(t, w(:,3), 'g');
hold off;
grid on;
legend('x', 'y', 'z');
title('Angular velocities');
xlabel('time [s]'); 
ylabel('angular rate [deg/s]');

set(fig2, 'Units', 'Inches');
pos1 = get(fig2, 'Position');
set(fig2, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos1(3), pos1(4)]);
print(fig2, '1_5_angular_velocities', '-dpdf', '-r0');

fig3 = figure (3); clf;
hold on;
plot(t, tau(:,1), 'b');
plot(t, tau(:,2), 'r');
plot(t, tau(:,3), 'g');
hold off;
grid on;
legend('x', 'y', 'z');
title('Control input');
xlabel('time [s]'); 
ylabel('input [Nm]');

set(fig3, 'Units', 'Inches');
pos1 = get(fig3, 'Position');
set(fig3, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos1(3), pos1(4)]);
print(fig3, '1_5_control_input', '-dpdf', '-r0');

fig4 = figure (4); clf;
hold on;
plot(t, phi_squiggle, 'b');
plot(t, theta_squiggle, 'r');
plot(t, psi_squiggle, 'g');

hold off;
grid on;
leg4 = legend('$\tilde{\phi}$', '$\tilde{\theta}$', '$\tilde{\psi}$');
set(leg4, 'Interpreter', 'latex');
title('Euler angle errors');
xlabel('time [s]');
ylabel('angle [deg]');

set(fig4, 'Units', 'Inches');
pos1 = get(fig4, 'Position');
set(fig4, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos1(3), pos1(4)]);
print(fig4, '1_5_euler_angle_errors', '-dpdf', '-r0');
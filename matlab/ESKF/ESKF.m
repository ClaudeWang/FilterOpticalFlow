function [state, test] = ESKF(sensor, varargin)
%ESTIMATE_VEL 6DOF velocity estimator
%   sensor - struct stored in provided dataset, fields include
%          - is_ready: logical, indicates whether sensor data is valid
%          - t: timestamp
%          - rpy, omg, acc: imu readings, you should not use these in this phase
%          - img: uint8, 240x376 grayscale image
%          - id: 1xn ids of detected tags
%          - p0, p1, p2, p3, p4: 2xn pixel position of center and
%                                four corners of detected tags
%            Y
%            ^ P3 == P2
%            | || P0 ||
%            | P4 == P1
%            o---------> X
%   varargin - any variables you wish to pass into the function, could be
%              a data structure to represent the map or camera parameters,
%              your decision. But for the purpose of testing, since we don't
%              know what inputs you will use, you have to specify them in
%              init_script by doing
%              estimate_vel_handle = ...
%                  @(sensor) estimate_vel(sensor, your personal input arguments);
%   vel - 3x1 velocity of the quadrotor in world frame
%   omg - 3x1 angular velocity of the quadrotor
persistent p q v b_a b_w err_p err_theta err_v err_ba err_bw t_prev t_interval Q sigma W R_noise

if isempty(p)
    
    % intialization.
    
    % noise
    % process noise.
    Q = diag([0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1])/5;
    
    % measurement noise.
    R_noise = 1e-10 * eye(3);
    
    % error state covariance
    sigma = diag([0.01 0.01 0.01 0.01 0.01 0.01 1e-3 1e-3 1e-3 ...
            0.025 0.025 0.005 0.001 0.001 0.001]);
       
    % nominal
    p = [0; 0; 0];
    v = [0; 0; 0];
    b_a = [0; 0; 0];
    b_w = [0; 0; 0];
    
    % error state
    err_p = [0; 0; 0];
    err_theta = [0; 0; 0];
    err_v = [0; 0; 0];
    err_ba = [0; 0; 0];
    err_bw = [0; 0; 0];
    
    t_prev = sensor.t;
    t_interval = 0.024;
    [~, ~] = estimate_vel_research(sensor);
    
    % solve for initial pose.
    % calculate rotation.
    acc = sensor.acc;
    acc_unit = acc / norm(acc);

    % This is used to formulate the first frame.
    % R_b_w = [r1, r2, r3].
    % align acc_unit and actual_g, the projection of z-axis of world is the
    % the third column of R_b_w.
    r3 = -acc_unit;

    % form r2 perpendicular to the r3.
    vector_random = [0; 0; 1];
    r2 = cross(r3, vector_random);
    r2 = r2 / norm(r2);

    % get r1 by crossing r2 and r3
    r1 = cross(r2, r3);
    R_b_w = [r1, r2, r3];
    
    R = R_b_w';
%     q = rotm2quat(R)';
    
    state = [p; v; q; b_a; b_w];
    
%     test = eul2rotm(sensor.rpy', 'XYZ');
    test = state;
    q = rotm2quat(eul2rotm(sensor.rpy', 'XYZ'))';
    return;
end

[vel, omg, R_vision] = estimate_vel_research(sensor);

% nominal state
t_interval = sensor.t - t_prev;
t_prev = sensor.t;

% state:
% [x, y, z, qw, qx, qy, qz, b_acc, b_omg]

% update nominal position
p = p + v * t_interval;

% update rotation.
omega = sensor.omg - b_w;

% conver the angular velocity to global.
w_norm = norm(omega);
angle = t_interval * w_norm;
axis = omega / w_norm;

delta_q = zeros(4, 1);
delta_q(1) = cos(angle / 2);
delta_q(2: 4) = sin(angle / 2) * axis;

q = quatmult(q, delta_q);
R = quat2rotm(q');

% R = eul2rotm(sensor.rpy', 'XYZ');
% test = R;

% update velocity.
v = v + (R * (sensor.acc - b_a) + [0; 0; -9.8]) * t_interval;

%% Error state update mean and covaraince.

% construct A matrix
A = zeros(15);
A(1:3, 4:6) = eye(3);
A(4:6, 7:9) = -R * skew(sensor.acc - b_a);
A(4:6, 10:12) = -R;
A(7:9, 13:15) = -R;

A = A * t_interval + eye(15);

% construct U matrix
U = zeros([15, 12]);

U(4:6, 1:3) = -R * t_interval^2;
U(7:9, 4:6) = -R * t_interval^2;
U(10: 15, 7:12) = eye(6) * t_interval;

% predict
old_err_state = [err_p; err_v; err_theta; err_ba; err_bw];
new_err_state =  A * old_err_state;
sigma = A * sigma * A' + U * Q * U';

%% Error state measurement model.
% measurement update, only use linear velocity as measurement. 3 * 15
% C = [zeros(3, 3), eye(3), zeros(3, 9)];

% H matrix.
rotation_q = rotm2quat([sqrt(2)/2, -sqrt(2)/2, 0; -sqrt(2)/2, -sqrt(2)/2, 0; 0, 0, -1] * quat2rotm(q'))';
v_quat = [0; v];
v_rotated = left_quat(quatmult(rotation_q, v_quat)) * diag([1, -1, -1, -1]) + right_quat(quatmult(v_quat, diag([1, -1, -1, -1]) * rotation_q));

H1 = [zeros(3), v_rotated(2:4, :), zeros(3), zeros(3), zeros(3)];
H2 = [eye(3), zeros(3, 12);
      zeros(3), eye(3), zeros(3, 9);
      zeros(4, 6), 0.5 * [-q(2), -q(3), -q(4); q(1), q(4), -q(3); -q(4), q(1), q(2); q(3), -q(2), q(1)], zeros(4, 6);
      zeros(3, 9), eye(3), zeros(3, 3);
      zeros(3, 12), eye(3)];
C = H1 * H2;
   
% kalman gain
K_p = (sigma * C') / (C * sigma * C' + R_noise);

% update the mean of error state
R_c_I = [sqrt(2)/2, -sqrt(2)/2, 0; -sqrt(2)/2, -sqrt(2)/2, 0; 0, 0, -1];
vel_in_world = R_c_I * quat2rotm(q') * v - R_c_I * skew(-R_c_I' * [-0.04; 0; -0.03]) * omega;
new_err_state = K_p * (vel - vel_in_world);

% update the covariance of the error state
sigma = (eye(15) - K_p * C) * sigma * (eye(15) - K_p * C)' ...
    + K_p * R_noise * K_p';

% update the state mean, give the final prediction.
p_final = p + new_err_state(1:3);
v_final = v + new_err_state(4:6);

dq = 0.5 * new_err_state(7:9);

if norm(dq) > 1
    dq = [1; dq] /  sqrt(1 + dq' * dq);
else
    dq = [sqrt(1 - dq' * dq); dq];
end

q_final = quatmult(dq, q);
b_a_final = b_a + new_err_state(10:12);
b_w_final = b_w + new_err_state(13:15);

% reset error state.
% G = blkdiag(eye(6), eye(3) + skew(0.5 * err_theta), eye(6));
G = eye(15);

err_p = [0; 0; 0];
err_theta = [0; 0; 0];
err_v = [0; 0; 0];
err_ba = [0; 0; 0];
err_bw = [0; 0; 0];

sigma = G * sigma * G';


state = [p_final; q_final; v_final; b_a_final; b_w_final];
% state = [p; q; v; b_a; b_w];
test = [p; q; v; b_a; b_w];

end

function quat = quatmult(q, r)

w1 = q(1);
x1 = q(2);
y1 = q(3);
z1 = q(4);

w2 = r(1);
x2 = r(2);
y2 = r(3);
z2 = r(4);

quat = [w1*w2-x1*x2-y1*y2-z1*z2;
        w1*x2+x1*w2+y1*z2-z1*y2;
        w1*y2-x1*z2+y1*w2+z1*x2;
        w1*z2+x1*y2-y1*x2+z1*w2
        ];

end

function T = skew(t)

T = [0, -t(1), t(2);
    t(1), 0, -t(3);
    -t(2), t(3), 0];

end

function Q = left_quat(q)
q0 = q(1);
q1 = q(2);
q2 = q(3);
q3 = q(4);

Q = [q0, -q1, -q2, -q3;
    q1, q0, -q3, q2;
    q2, q3, q0, -q1;
    q3, -q2, q1, q0];
end

function Q = right_quat(q)
q0 = q(1);
q1 = q(2);
q2 = q(3);
q3 = q(4);

Q = [q0, -q1, -q2, -q3;
    q1, q0, q3, -q2;
    q2, -q3, q0, -q1;
    q3, q2, -q1, q0;];
end




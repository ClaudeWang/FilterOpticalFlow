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
persistent p q v b_a b_w err_p err_q err_v err_ba err_bw t_prev t_interval

if isempty(p)
    
    % intialization.
    
    % nominal
    p = [0; 0; 0];
    v = [0; 0; 0];
    b_a = [0; 0; 0];
    b_w = [0; 0; 0];
    
    % error state
    err_p = [0; 0; 0];
    err_q = [0; 0; 0; 0];
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
    
    state = [p; q; b_a; b_w];
    
    test = eul2rotm(sensor.rpy', 'XYZ');
    q = rotm2quat(test)';
    return;
end

[vel, omg, R_vision] = estimate_vel_research(sensor);

% if size(vel, 1) == 0
%     return;
% end

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
% 
delta_q = zeros(4, 1);
delta_q(1) = cos(angle / 2);
delta_q(2: 4) = sin(angle / 2) * axis;
% 
q = quatmult(delta_q, q);
R = quat2rotm(q');

% q = rotm2quat(R)';


% update velocity.
% sensor.acc;
% quat2rotm(q') * (sensor.acc - b_a);
% v = v + (R * (sensor.acc - b_a) + [0; 0; 9.8]) * t_interval;

state = [p; q; v; b_a; b_w];
test = R;
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






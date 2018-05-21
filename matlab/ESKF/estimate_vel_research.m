function [vel, omg, R_test] = estimate_vel_research(sensor, varargin)
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
persistent point_tracker last_points last_t last_vel last_omg last_test_R

if size(sensor.id, 2) == 0
    vel = [];
    omg = [];
    flow = [];
    depth = [];
    points = [];
    return;
end

if isempty(last_vel)
    last_vel = [0; 0; 0];
    last_omg = [0; 0; 0];
    last_test_R = eye(3);
end


if isempty(point_tracker)
    minQuality = 0.35;
    nPyrLevels = 1;
    maxBiError = 3;
    
    img0 = sensor.img;
    new_points_feature = detectMinEigenFeatures(img0,'MinQuality',minQuality);
    new_points_feature = new_points_feature.Location;

    if size(new_points_feature, 1) > 500
        new_points_feature = new_points_feature(1:500, :);
    end
        
    point_tracker = vision.PointTracker('NumPyramidLevels',nPyrLevels,...
                                       'MaxBidirectionalError',maxBiError);
    initialize(point_tracker, new_points_feature, img0);
    last_points = new_points_feature;
    last_t = sensor.t;
    vel = last_vel;
    omg = last_omg;
    flow = [];
    points = [];
    depth = [];
    R_test = last_test_R;
    return
end

% all the y and x indices
id = sensor.id;
y_index = floor(id / 12);
x_index = mod(id, 12);
x_left_up = 0.152 * 2 * x_index;
y_left_up = 0.152 * 2 * y_index;

y_left_up(y_index >= 3) = y_left_up(y_index >=3) + 0.026;
y_left_up(y_index >= 6) = y_left_up(y_index >= 6) + 0.026;

x_center = x_left_up + 0.152 / 2;
y_center = y_left_up + 0.152 / 2;

x_right_down = x_left_up + 0.152;
y_right_down = y_left_up + 0.152;

x_right_up = x_left_up;
y_right_up = y_left_up + 0.152;

x_left_down = x_left_up + 0.152;
y_left_down = y_left_up;


% put all the x together.
x_world = [x_left_up; x_center; x_right_down; x_right_up; x_left_down];
y_world = [y_left_up; y_center; y_right_down; y_right_up; y_left_down];

% Initialize pixels.
all_sensor_x = [sensor.p4(1, :);
    sensor.p0(1, :);
    sensor.p2(1, :);
    sensor.p3(1, :);
    sensor.p1(1, :)];

all_sensor_y = [sensor.p4(2, :);
    sensor.p0(2, :);
    sensor.p2(2, :);
    sensor.p3(2, :);
    sensor.p1(2, :)];

april_tag_x = x_world(:);
april_tag_y = y_world(:);

pixel_x = all_sensor_x(:);
pixel_y = all_sensor_y(:);

zero_entries = zeros(5 * size(sensor.id, 2), 1);
one_entries = ones(5 * size(sensor.id, 2), 1);

Ax = [-april_tag_x, zero_entries, pixel_x.*april_tag_x, -april_tag_y, zero_entries, pixel_x.*april_tag_y,... 
    -one_entries, zero_entries, pixel_x];
Ay = [zero_entries, -april_tag_x, pixel_y.*april_tag_x, zero_entries, -april_tag_y, pixel_y.*april_tag_y,...
    zero_entries, -one_entries, pixel_y];


A = [Ax; Ay];

% perform svd to estimate the best projection homography.
[~,~,V] = svd(A);
h = V(:, end);
h = h / h(9);
H = reshape(h, [3, 3]);

% intrinsic matrix to calculate the pose.
K = [311.0520 0        201.8724;
    0         311.3885 113.6210;
    0         0        1];
pose_c_w = K \ H;

% convert the orientation from "camera to world" to "world to camera".
r1 = pose_c_w(:, 1);
r2 = pose_c_w(:, 2);
r3 = cross(r1, r2);
t = 2 * pose_c_w(:, 3) / (norm(r1) + norm(r2));

R = [r1, r2, r3];
[U, ~, V] = svd(R);
R = U * [1, 0, 0; 0, 1, 0; 0, 0, det(U*V')] * V';
% the rotation of the camera.
R_c_w = R;
t_c_w = t;

R_c_I = [sqrt(2)/2, -sqrt(2)/2, 0;
        -sqrt(2)/2, -sqrt(2)/2, 0;
        0, 0, -1];
R_test = R_c_w' * R_c_I;
last_test_R = R_test;

% track points
new_image = sensor.img;
[new_points,valid] = step(point_tracker, new_image);
num_points_left = nnz(valid);

% all the points can be converted to the world frame.
inv_H = inv(H);
inv_H = inv_H / inv_H(3, 3);
old_points_feature = last_points(valid, :);
new_points_feature = new_points(valid, :);
last_points = new_points;

time_interval = sensor.t - last_t;
last_t = sensor.t;

new_points_feature_world = inv_H * [new_points_feature'; ones(1, size(new_points_feature, 1))];
new_points_feature_world = bsxfun(@rdivide,new_points_feature_world, new_points_feature_world(3, :));

% now all the points can be converted back to the camera frame.
new_points_feature_world(3, :) = zeros([1, num_points_left]);
new_points_feature_camera = R_c_w * new_points_feature_world + t_c_w;
Z = new_points_feature_camera(3, :);
Z = Z';

% figure(1)
% imshow(sensor.img)
% hold on;
% scatter(new_points_feature(:, 1), new_points_feature(:, 2))
% drawnow;
% hold off;


ideal_coordinates_new = inv(K) * [new_points_feature'; ones(1, num_points_left)];
x = ideal_coordinates_new(1, :)';
y = ideal_coordinates_new(2, :)';
ideal_coordinates_old = inv(K) * [old_points_feature'; ones(1, num_points_left)];
% optical_flow = (ideal_coordinates_new - ideal_coordinates_old) / time_interval;
optical_flow = (ideal_coordinates_new - ideal_coordinates_old) / 0.0205;

flow = optical_flow(1:2, :);
points = new_points_feature(:, 1:2)';
depth = Z';

zeros_entries = zeros([num_points_left, 1]);
all_A = [-1./Z, zeros_entries, x./Z, x.*y, -(1 + x.^2), y;
     zeros_entries, -1./Z, y./Z, 1 + y.^2, -x.*y, -x];
% 
all_b = [optical_flow(1, :)'; optical_flow(2, :)'];
% 
% result = A \ b;

% RANSAC to reject outliers.
tolerance = 0.3;
max_inlier = -1;
max_inlier_mask = [];
for i = 1: 100
    % randomly choose three points.
    index = randsample(num_points_left, 3);
    sample_x = x(index);
    sample_y = y(index);
    sample_Z = Z(index);
    zeros_entries = zeros([3, 1]);
    A = [-1./sample_Z, zeros_entries, sample_x./sample_Z, sample_x.*sample_y, -(1 + sample_x.^2), sample_y;
        zeros_entries, -1./sample_Z, sample_y./sample_Z, 1 + sample_y.^2, -sample_x.*sample_y, -sample_x];
    b = [optical_flow(1, index)'; optical_flow(2, index)'];
    result = A \ b;
    
    diff_opt = all_A * result - all_b;
    x_diff = diff_opt(1: num_points_left);
    y_diff = diff_opt(num_points_left + 1: end);
    
    mask = sqrt(x_diff.^2 + y_diff.^2) < tolerance;
    inlier_num = nnz(mask);
    if inlier_num > max_inlier
        max_inlier  = inlier_num;
        max_inlier_mask = mask;
%         if inlier_num > num_points_left * 0.95
%             break
%         end
    end
end

inlier_x = x(max_inlier_mask);
inlier_y = y(max_inlier_mask);
inlier_z = Z(max_inlier_mask);
zeros_entries = zeros([max_inlier, 1]);

A = [-1./inlier_z, zeros_entries, inlier_x./inlier_z, inlier_x.*inlier_y, -(1 + inlier_x.^2), inlier_y;
    zeros_entries, -1./inlier_z, inlier_y./inlier_z, 1 + inlier_y.^2, -inlier_x.*inlier_y, -inlier_x];
b = [optical_flow(1, max_inlier_mask)'; optical_flow(2, max_inlier_mask)'];
best_result = A \ b;



% if size(result, 1) == 0
%     best_result = all_A \ all_b;
% end


% reinitialize if the number of points left is too low.
left_points_threshold = 300;
if num_points_left < left_points_threshold
    clear point_tracker
end

vel = best_result(1: 3);
omg = best_result(4: 6);

omg = R_c_w' * omg;
vel = R_c_w' * vel;
if ~isempty(last_vel)
    vel = 0.5 * vel + 0.5 * last_vel;
end
last_vel = vel;
last_omg = omg;
end

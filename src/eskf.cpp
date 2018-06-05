//
// Created by Haumin Zhu on 5/25/18.
//

#include "include/eskf/eskf.h"

Matrix4f leftquat(Vector4f q);
Matrix4f rightquat(Vector4f q);

void ESKF::filter(float dt, MatrixXf imu, MatrixXf flow, MatrixXf K, MatrixXf depth, MatrixXf position, bool is_ready){
    prediction(dt, imu);
    if(is_ready){
        measurement(flow, K, depth, position, imu);
    }
}


void ESKF::prediction(float dt, MatrixXf imu){
    std::cout << "dt: "<< dt << std::endl;
    std::cout << "imu: "<< imu << std::endl;

    Vector3f p = _nom_state.head(3);
    Vector3f v = _nom_state.segment(3, 3);
    Vector4f q = _nom_state.segment(6, 4);

    Vector3f acc = imu.row(0).head(3);
    Vector3f omg = imu.row(0).tail(3);

    p = p + v * dt;

    Vector3f axis = omg.normalized();
    float angle = omg.norm() * dt;

    Vector4f dq;
    Vector3f second = sin(angle / 2) * axis;
    dq = {(float)(cos(angle / 2)), second[0], second[1], second[2]};

    q = quatMultiply(q, dq);

    Quaternionf quat(q[0], q[1], q[2], q[3]);

    Matrix3f wRb = quat.normalized().toRotationMatrix();

    Vector3f g = {0, 0, float(-9.8)};
    Vector3f bias_acc = _nom_state.segment(10, 3);

    Vector3f tmp = wRb * (acc - bias_acc);
    v = v + (tmp + g) * dt;

    MatrixXf A = MatrixXf::Zero(15, 15);
    A.block(0, 3, 3, 3) = MatrixXf::Identity(3, 3);
    A.block(3, 6, 3, 3) = -wRb * skew(acc - _nom_state.segment(10, 3));

    A.block(3, 9, 3, 3) = -wRb;
    A.block(6, 12, 3, 3) = -wRb;


    A = A * dt + MatrixXf::Identity(15, 15);

    MatrixXf U = MatrixXf::Zero(15, 12);
    U.block(3, 0, 3, 3) = -wRb * dt * dt;
    U.block(6, 3, 3, 3) = -wRb * dt * dt;
    U.block(9, 6, 6, 6) = MatrixXf::Identity(6, 6) * dt;

    _err_state = A * _err_state;
    _cov = A * _cov * A.transpose() + U * _Q * U.transpose();

    VectorXf new_nom(16);
    new_nom << p[0], p[1], p[2], v[0], v[1], v[2], q[0], q[1], q[2], q[3], _nom_state[10], _nom_state[11],
            _nom_state[12], _nom_state[13], _nom_state[14], _nom_state[15];
    _nom_state = new_nom;

    return;
}

void ESKF::measurement(MatrixXf flow, MatrixXf K, MatrixXf depth, MatrixXf position, MatrixXf imu){

    VectorXf vel = estimateVelocity_.velocityFromRANSAC(_ransac_iteration, _ransac_threshold, flow, K, depth, position).col(0).head(3);

    Vector4f q = _nom_state.segment(6, 4);
    Quaternionf quat(q[0], q[1], q[2], q[3]);
    Matrix3f wRb = quat.normalized().toRotationMatrix();

    Matrix3f rot2q;
    rot2q << sqrt(2)/2, -sqrt(2)/2, 0, -sqrt(2)/2, -sqrt(2)/2, 0, 0, 0, -1;
    rot2q = (rot2q * wRb).transpose();
    Quaternionf rot_q(rot2q);
    Vector4f rotation_q;

    rotation_q << rot_q.w(), rot_q.x(), rot_q.y(), rot_q.z();

    Vector4f v_quat;
    Vector3f v = _nom_state.segment(3, 3);
    v_quat << 0, v[0], v[1], v[2];

    Vector4f dia_tmp;
    dia_tmp << 1, -1, -1, -1;
    Matrix4f v_rotated;
    v_rotated = leftquat(quatMultiply(rotation_q, v_quat)) * dia_tmp.asDiagonal() + rightquat(quatMultiply(v_quat, dia_tmp.asDiagonal() * rotation_q));

    Matrix3f R_2_1;
    R_2_1 << sqrt(2)/2, -sqrt(2)/2, 0, -sqrt(2)/2, -sqrt(2)/2, 0, 0, 0, -1;
    Matrix3f R_1_0 = wRb.transpose();

    Vector3f d;
    d << -0.04, 0, -0.03;
    Matrix3f j_w_b = R_2_1 * skew(-R_2_1.transpose() * d);


    MatrixXf H1 = MatrixXf::Zero(3, 16);
    MatrixXf H2 = MatrixXf::Zero(16, 15);

    H1.block(0, 3, 3, 3) = R_2_1 * R_1_0;
    H1.block(0, 6, 3, 4) = v_rotated.block(1, 0, 3, 4);
    H1.block(0, 13, 3, 3) = j_w_b;

    H2.block(0, 0, 6, 6) = MatrixXf::Identity(6, 6);
    MatrixXf q_op(4, 3);
    q_op << -q[1], -q[2], -q[3], q[0], q[3], -q[2], -q[3], q[0], q[1], q[2], -q[1], q[0];
    H2.block(6, 6, 4, 3) = 0.5 * q_op;
    H2.block(10, 9, 6, 6) = MatrixXf::Identity(6, 6);

    MatrixXf C(3, 15);
    C = H1 * H2;

    MatrixXf Kp(15, 3);
    Kp = _cov * C.transpose() * (C * _cov * C.transpose() + _R).inverse();

    std::cout << Kp << std::endl;

    Vector3f omg = imu.row(0).tail(3);
    std::cout << omg << std::endl;

    Vector3f vel_in_cam = R_2_1 * wRb * v - R_2_1 * skew(-R_2_1.transpose() * d) * omg;
    std::cout << vel_in_cam << std::endl;

    _err_state = Kp * (vel - vel_in_cam);

    std::cout << _err_state << std::endl;

    _cov = (MatrixXf::Identity(15, 15) - Kp * C) * _cov * (MatrixXf::Identity(15, 15) - Kp*C).transpose() + Kp * _R * Kp.transpose();

    Vector3f p_true = _nom_state.head(3) + _err_state.head(3);
    Vector3f v_true = _nom_state.segment(3, 3) + _err_state.segment(3, 3);

    Vector3f d_omg = 0.5 * _err_state.segment(6, 3);

    std::cout << d_omg << std::endl;
    Vector4f dq;

    if(d_omg.norm() > 1){
        Vector4f tmp;
        tmp << 1, d_omg(0), d_omg(1), d_omg(2);
        dq = tmp / (float)sqrt(1 + d_omg.dot(d_omg));
    } else {
        dq = {(float)sqrt(1 - d_omg.dot(d_omg)), d_omg[0], d_omg[1], d_omg[2]};
    }
    std::cout << dq << std::endl;

    Vector4f q_true = quatMultiply(dq, _nom_state.segment(6, 4));
    Vector3f b_a_true = _nom_state.segment(10, 3) + _err_state.segment(9, 3);
    Vector3f b_w_true = _nom_state.tail(3) + _err_state.tail(3);

    MatrixXf G = MatrixXf::Identity(15, 15);
    _err_state.setZero();

    _cov = G * _cov * G.transpose();

    _true_state << p_true(0), p_true[1], p_true[2], v_true[0], v_true[1], v_true[2],
            q_true[0], q_true[1], q_true[2], q_true[3], b_a_true[0], b_a_true[1], b_a_true[2], b_w_true[0], b_w_true[1], b_w_true[2];

    std::cout << _true_state << std::endl;
    return;
}

Vector4f ESKF::quatMultiply(Vector4f q, Vector4f r){
    Vector4f result(q[0] * r[0] - q[1] * r[1] - q[2] * r[2] - q[3] * r[3],
                    q[0] * r[1] + q[1] * r[0] + q[2] * r[3] - q[3] * r[2],
                    q[0] * r[2] - q[1] * r[3] + q[2] * r[0] + q[3] * r[1],
                    q[0] * r[3] + q[1] * r[2] - q[2] * r[1] + q[3] * r[0]);
    return result;
}

Matrix3f ESKF::skew(Vector3f vec){
    Matrix3f s;
    s << 0, -vec[0], vec[1], vec[0], 0, -vec[2], -vec[1], vec[2], 0;
    return s;
}

Matrix4f leftquat(Vector4f q){
    Matrix4f Q;
    Q << q[0], -q[1], -q[2], -q[3], q[1], q[0], -q[3], q[2], q[2], q[3], q[0], -q[1], q[3], -q[2], q[1], q[0];
    return Q;
}

Matrix4f rightquat(Vector4f q){
    Matrix4f Q;
    Q << q[0], -q[1], -q[2], -q[3], q[1], q[0], q[3], -q[2], q[2], -q[3], q[0], -q[1], q[3], q[2], -q[1], q[0];
    return Q;
}
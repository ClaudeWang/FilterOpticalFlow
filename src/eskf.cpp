//
// Created by Haumin Zhu on 5/25/18.
//

#include "include/eskf/eskf.h"


void ESKF::prediction(float dt, VectorXf imu){
    Vector3f p = _nom_state.head(3);
    Vector4f q = _nom_state.segment(3, 7);
    Vector3f v = _nom_state.segment(7, 10);

    Vector3f omg = imu.head(3);
    Vector3f acc = imu.tail(3);

    p = p + v * dt;

    Vector3f axis = omg.normalized();
    float angle = omg.norm() * dt;

    Vector4f dq;
//    dq【0] = cos(angle / 2);
    Vector3f second = sin(angle / 2) * axis;
    dq = {(float)(cos(angle / 2)), second[0], second[1], second[2]};
//    dq.tail(3) = sin(angle / 2) * axis;

    q = quatMultiply(q, dq);

    Quaternionf quat(q[0], q[1], q[2], q[3]);

    Matrix3f wRb = quat.normalized().toRotationMatrix();

    Vector3f g = {0, 0, float(-9.8)};
    Vector3f bias_acc = _nom_state.segment(10, 13);

    Vector3f tmp = wRb * (acc - bias_acc);
    v = v + (tmp + g) * dt;

//    Matrix<float, 15, 15> A;
    MatrixXf A(15, 15);
    A.block(0, 3, 3, 3) = MatrixXf::Identity(3, 3);
    A.block(3, 6, 3, 3) = -wRb * skew(acc - _nom_state.segment(10, 13));
    A.block(3, 9, 3, 3) = -wRb;
    A.block(6, 12, 3, 3) = -wRb;

    A = A * dt + MatrixXf::Identity(15, 15);

//    Matrix<float, 15, 12> U;
    MatrixXf U(15, 12);
    U.block(3, 0, 3, 3) = -wRb * dt * dt;
    U.block(6, 3, 3, 3) = -wRb * dt * dt;
    U.block(9, 6, 6, 6) = MatrixXf::Identity(6, 6) * dt;

    _err_state = A * _err_state;
    _cov = A * _cov * A.transpose() + U * _Q * U.transpose();

    return;
}

void ESKF::measurement(MatrixXf flow, MatrixXf K, MatrixXf depth, MatrixXf position){

    VectorXf vel = estimateVelocity_.velocityFromRANSAC(_ransac_iteration, _ransac_threshold, flow, K, depth, position).col(0).head(3);

//    Matrix<float, 3, 15> C;
    MatrixXf C(3, 15);
    C.block(0, 3, 3, 3) = MatrixXf::Identity(3, 3);

//    Matrix<float, 15, 3> Kp;
    MatrixXf Kp(15, 3);
    Kp = _cov * C.transpose() * (C * _cov * C.transpose() + _R).inverse();
    _err_state = Kp * (vel - _nom_state.segment(7, 10));

    _cov = (MatrixXf::Identity(15, 15) - Kp * C) * _cov * (MatrixXf::Identity(15, 15) - Kp*C).transpose() + Kp * _R * Kp.transpose();

    Vector3f p_true = _nom_state.head(3) + _err_state.head(3);
    Vector3f v_true = _nom_state.segment(7, 10) + _err_state.segment(3, 6);

    Vector3f d_omg = 0.5 * _err_state.segment(6, 9);
    Vector4f dq;

    if(d_omg.norm() > 1){
        Vector4f tmp;
        tmp << 1, d_omg(0), d_omg(1), d_omg(2);
        dq = tmp / (float)sqrt(1 + d_omg.dot(d_omg));
    } else {
        dq = {(float)sqrt(1 - d_omg.dot(d_omg)), d_omg[0], d_omg[1], d_omg[2]};
    }

    Vector4f q_true = quatMultiply(dq, _nom_state.segment(3, 7));
    Vector3f b_a_true = _nom_state.segment(10, 13) + _err_state.segment(9, 12);
    Vector3f b_w_true = _nom_state.tail(3) + _err_state.tail(3);

    MatrixXf G = MatrixXf::Identity(15, 15);
    _err_state.setZero();

    _cov = G * _cov * G.transpose();

//    _true_state(0) = p_true(0);
    _true_state << p_true(0), p_true[1], p_true[2], q_true[0], q_true[1], q_true[2], q_true[3], v_true[0], v_true[1],
            v_true[2], b_a_true[0], b_a_true[1], b_a_true[2], b_w_true[0], b_w_true[1], b_w_true[2];
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
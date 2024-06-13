#ifndef FASTER_LIO_INVKF_H
#define FASTER_LIO_INVKF_H
#include <Eigen/Dense>


#include "Eigen/src/Core/Matrix.h"
#include "Eigen/src/Geometry/Quaternion.h"
#include "liepp/SEn3.h"
#include "liepp/SO3.h"


namespace InvariantKF{
    template <typename T>
    inline Eigen::Matrix<T, 3, 3> SKEW_SYM_MATRIX(const Eigen::Matrix<T, 3, 1> &v) {
        Eigen::Matrix<T, 3, 3> m;
        m << 0.0, -v[2], v[1], v[2], 0.0, -v[0], -v[1], v[0], 0.0;
        return m;
    } 
    template <typename T>
    Eigen::Matrix<T, 3, 3> Exp(const Eigen::Matrix<T, 3, 1> &&ang) {
        T ang_norm = ang.norm();
        Eigen::Matrix<T, 3, 3> Eye3 = Eigen::Matrix<T, 3, 3>::Identity();
        if (ang_norm > 0.0000001) {
            Eigen::Matrix<T, 3, 1> r_axis = ang / ang_norm;
            Eigen::Matrix<T, 3, 3> K;
            K = SKEW_SYM_MATRIX(r_axis);
        /// Roderigous Tranformation
            return Eye3 + std::sin(ang_norm) * K + (1.0 - std::cos(ang_norm)) * K * K;
        } else {
            return Eye3;
        }
    }  
    // template <typename state, int process_noise_dof, typename input, typename measurement ,
    //       int measurement_noise_dof = 0>
    class invkf{
        public:
        invkf(){
            
            // imu_state.setIdentity();
            // gravity.setZero();
            // theta.setZero();
            // exrinsic.setIdentity();
            // P.setIdentity();
            // Q.setIdentity();
            
            // P(6, 6) = P(7, 7) =  P(8, 8) = 0.00001;
            // P(9, 9) =  P(10, 10) =  P(11, 11) = 0.00001;
            // P(15, 15) =  P(16, 16) =  P(17, 17) = 0.0001;
            // P(18, 18) =  P(19, 19) =  P(20, 20) = 0.001;
            // P(21, 21) =  P(22, 22) = 0.00001;
            
        };
        struct State24{
            liepp::SEn3<2> imu_state; //9
            Eigen::Vector3d gravity; //3
            //Eigen::VectorXd theta;
            Eigen::Vector3d bg;  //3
            Eigen::Vector3d ba;  //3
            Eigen::Quaterniond ext_r;
            Eigen::Vector3d ext_t;
            State24(){
                imu_state.setIdentity();//9
                gravity.setZero();//3
                bg.setZero();//3
                ba.setZero();//3
                ext_r.setIdentity();//3
                ext_t.setZero();//3
            }
        };
        void predict(const double dt, const Eigen::Vector3d &acc_, const Eigen::Vector3d &gyr_){
            acc = acc_- state.ba; 
            gyr = gyr_- state.bg;  
            
            auto rotation = (state.imu_state.R  * liepp::SO3<double>().exp(gyr * dt)); // r * exp(w * dt)
            //Eigen::Quaterniond rotation_qua = rotation.asQuaternion();
            Eigen::Vector3d velocity = state.imu_state.x[0] + (rotation * acc + state.gravity) * dt; // v+(r * a+ g) * dt
            Eigen::Vector3d position = state.imu_state.x[1] * dt + state.imu_state.x[1] + 0.5 * (rotation * acc + state.gravity) * dt * dt; //p + v * dt + 0.5 *()
            Eigen::Matrix3d rotation_mat = rotation.asMatrix();
            std::array<Eigen::Vector3d, 2> temp{velocity,position};
            liepp::SEn3<2> new_state(rotation,temp);//set state
            state.imu_state = new_state;

            Eigen::MatrixXd A = Eigen::MatrixXd::Zero(24,24);
            A.block(3,0,3,3) = SKEW_SYM_MATRIX(state.gravity);
            A.block(6,3,3,3) = Eigen::Matrix3d::Identity();
            A.block(0,9,3,3) = rotation_mat * (-1);
            A.block(3,9,3,3) = -1 * SKEW_SYM_MATRIX(velocity) * rotation_mat;
            A.block(3,12,3,3) = -1 * rotation_mat;
            A.block(6,9,3,3) = -1 * SKEW_SYM_MATRIX(position) * rotation_mat;

            //discretization
            
            Eigen::MatrixXd I = Eigen::MatrixXd::Identity(24,24);
            //approximation of exp (A * dt) , msckf , need to test
            Eigen::MatrixXd Adt = A * dt;
            Eigen::MatrixXd Adt_square = Adt * Adt;
            Eigen::MatrixXd Adt_cube = Adt_square * Adt;
            Eigen::MatrixXd phi = I + Adt + 0.5 * Adt_square + (1.0/6.0) * Adt_cube;
            Eigen::MatrixXd Adj = I;
            Adj.block(0,0,9, 9) = state.imu_state.Adjoint();
            Eigen::MatrixXd PhiAdj = phi * Adj;
            Eigen::MatrixXd Qk_hat = PhiAdj * Q * PhiAdj.transpose() * dt;
            P =  phi * P *phi.transpose() + Qk_hat;
            
            
            return;
        }
        void update(){
            liepp::SEn3<2> state_ = state.imu_state;
            Eigen::MatrixXd K;
            Eigen::MatrixXd H;
            Eigen::MatrixXd z;
            calcZH(state_, z,H);
        }
        void calcZH(liepp::SEn3<2> &state, Eigen::MatrixXd &z, Eigen::MatrixXd &H ){

        }

        void setP(Eigen::Matrix<double, 24, 24> P_){
            P = P_;
        }
        const Eigen::Matrix<double, 24, 24>& getP(){
            return P;
        }
        void setQ(Eigen::Matrix<double, 24, 24> Q_){
            Q = Q_;
        }
        Eigen::Vector3d getGravity(){
            return state.gravity;
        }
        void setX(const State24 &state_){
            state = state_;
        }
        const State24& getX(){
            return state;
        };

        private:
        Eigen::Vector3d acc,gyr;
        // liepp::SEn3<2> imu_state;  //rvp
        // Eigen::Vector3d gravity;   //gravity
        // Eigen::VectorXd theta;     //bias
        // Eigen::Quaterniond ext_r;
        // Eigen::Vector3d ext_t;
        // liepp::SEn3<1> exrinsic; 
        State24 state;  
        Eigen::Matrix<double, 24, 24> P; 
        Eigen::Matrix<double, 24 ,24> Q;

        
        // state imu_state;
        
        // input imu_input;
        
    };
}
#endif
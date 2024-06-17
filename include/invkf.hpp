#ifndef FASTER_LIO_INVKF_H
#define FASTER_LIO_INVKF_H
#include <fcntl.h>
#include <Eigen/Dense>
#include <iostream>
#include <memory>
#include <fstream>

#include "Eigen/src/Core/Matrix.h"
#include "Eigen/src/Core/util/Constants.h"
#include "Eigen/src/Geometry/Quaternion.h"

#include "liepp/SEn3.h"
#include "liepp/SO3.h"
#include "use-ikfom.hpp"


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

    template <typename T>
    struct LioZHModel {
        bool valid;
        bool converge;
        Eigen::Matrix<T, Eigen::Dynamic, 1> z;
        Eigen::Matrix<T, Eigen::Dynamic, 1> h;
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> h_v;
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> h_x;
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> R;
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

        using measurementModel_dyn_share = std::function<void(State24 &, LioZHModel<double> &)>;
        measurementModel_dyn_share obsModel;
        void initObsModel(measurementModel_dyn_share obsModel_){
            obsModel = obsModel_;
        }
        void predict(const double dt, const Eigen::Vector3d &acc_, const Eigen::Vector3d &gyr_){
            acc = acc_- state.ba; 
            gyr = gyr_- state.bg;  
            
            auto rotation = ( state.imu_state.R  * faster_lio::SO3::exp(gyr * dt) ); // r * exp(w * dt)
            std::cout << state.gravity(0) <<" "<< state.gravity(1) <<" " << state.gravity(2) << std::endl<<std::endl;
            Eigen::Vector3d velocity = state.imu_state.x[0] + (rotation * acc + state.gravity) * dt; // v+(r * a+ g) * dt
            Eigen::Vector3d position = state.imu_state.x[0] * dt + state.imu_state.x[1] + 0.5 * (rotation * acc + state.gravity) * dt * dt; //p + v * dt + 0.5 *()
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
            //A.block(15,15,9,9) = Eigen::MatrixXd::Identity(9,9);

            //discretization
            
            Eigen::MatrixXd I = Eigen::MatrixXd::Identity(24,24);
            //approximation of exp (A * dt) , msckf , need to test
            Eigen::MatrixXd Adt = A * dt;
            Eigen::MatrixXd Adt_square = Adt * Adt;
            Eigen::MatrixXd Adt_cube = Adt_square * Adt;
            Eigen::MatrixXd phi = I + Adt + 0.5 * Adt_square + (1.0/6.0) * Adt_cube;
            Eigen::MatrixXd Adj = I;
            //std::cout <<" test"<< phi - Adt << std::endl;
            Adj.block(0,0,9, 9) = state.imu_state.Adjoint();
            Eigen::MatrixXd PhiAdj = phi * Adj;
            Eigen::MatrixXd Qk_hat = PhiAdj * Q * PhiAdj.transpose() * dt;
            P =  phi * P *phi.transpose() + Qk_hat;
            
            
            return;
        }
        void update(double R, double &Sovel_H_time){
            State24 state_prop = state;
            Eigen::MatrixXd P_prop = P;
            LioZHModel<double> liozh_model;
            liozh_model.valid = true;
            liozh_model.converge = true;
            int t = 0;
            Eigen::Matrix<double, 24, 1> K_h;
            Eigen::Matrix<double, 24, 24> K_x;
            Eigen::Matrix<double, 24, 1> dx_new =  Eigen::Matrix<double, 24, 1>::Zero();
            //先不做迭代
            liozh_model.valid = true;
            obsModel(state_prop,liozh_model);
            if(!liozh_model.valid){
                return;
            }
            Eigen::Matrix<double, Eigen::Dynamic, 24> h_x_ = liozh_model.h_x;
            Eigen::Matrix<double, 24, 1> dx = errState(state_prop);//误差
            //std::cout << "dx" << dx << std::endl;
            /**
             * @brief 卡尔曼： 
             K = P * HT (H * P * HT + R) ^-1
             K = (HT * R^-1 * H  + P^-1)^-1 * HT * R^-1
             x = exp(kz) x
             * K = (H_t*H+(P_/0.001).inverse()).inverse()*H_t;
             */
             
            //K_x = (h_x_.transpose() * h_x_ + (P/0.001).inverse()).inverse() * h_x_.transpose();
            Eigen::MatrixXd HTH = h_x_.transpose() * h_x_;
            Eigen::MatrixXd K;
            K = (HTH +  (P/0.001).inverse()).inverse() * h_x_.transpose();
            K_x = K * h_x_;

            // Eigen::MatrixXd P_temp = (P / 0.001).inverse(); 
            // Eigen::MatrixXd HTH = h_x_.transpose() * h_x_;
            // P_temp += HTH;//hx现在是完整的24维矩阵，不需要分块了
            // Eigen::MatrixXd  P_inv = P_temp.inverse();
            // K_h = P_inv * h_x_.transpose() * liozh_model.h; //kh
            // K_x.setZero();  // = cov_::Zero();
            // K_x = P_inv * HTH;
            // // Eigen::Matrix<double, 24, 1> dx_ =  K_h + (K_x - Eigen::Matrix<double, 24, 24>::Identity()) * dx ;
            // Eigen::Matrix<double, 24, 1> dx_ =  -1 *  (K_h + (K_x - Eigen::Matrix<double, 24, 24>::Identity()) * dx) ;//第一次dx_new为
            // std::cout << liozh_model.h.size() << std::endl;   
            Eigen::Matrix<double, 24, 1> dx_ =   K * liozh_model.h ;

            std::cout << "dx" << dx_ << std::endl;
            state.imu_state = liepp::SEn3<2>::exp(dx_.segment(0,9)) * state.imu_state ;

            
            state.bg = dx_.segment(9,3) + state.bg;
            state.ba = dx_.segment(12,3) + state.ba;
            state.gravity = dx_.segment(15,3) + state.gravity;
            state.ext_r = (liepp::SO3<double>::exp(dx_.segment(18,3)) * state.ext_r).asQuaternion();
            state.ext_t = dx_.segment(21,3) + state.ext_t;
            P = (Eigen::MatrixXd::Identity(24,24) - K_x ) * P;



        }                       
        const Eigen::Matrix<double,24,1> errState (State24 state_){
            Eigen::Matrix<double,24,1> result;

            result.segment(0,9) = liepp::SEn3<2>::log(state.imu_state * state_.imu_state.inverse() );
           
            result.segment(9,3) = state.bg - state_.bg;
            result.segment(12,3) = state.ba - state_.ba;
            result.segment(15,3) = state.gravity - state_.gravity;
            result.segment(18,3) = liepp::SO3<double>::log(state.ext_r * state_.ext_r.inverse());
            result.segment(21,3) = state.ext_t -state_.ext_t;
            return result;
        }
        // const State24 getErrorState(const State24 &state_hat, const State24 &state_gt){

        //     //err = state_hat - state_gt = 
        //     State24 err_state;
        //     err_state.imu_state = state_hat.imu_state * state_gt.imu_state.inverse();
        //     err_state.ba = state_hat.ba - state_gt.ba;
        //     err_state.bg = state_hat.bg - state_gt.bg;
        //     err_state.ext_r = state_hat.ext_r * state_gt.ext_r.inverse();
        //     err_state.ext_t = state_hat.ext_t - state_gt.ext_t;
        //     return err_state;
        // }
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
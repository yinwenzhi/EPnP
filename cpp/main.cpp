// Copyright (c) 2009, V. Lepetit, EPFL
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met: 

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer. 
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution. 

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// The views and conclusions contained in the software and documentation are those
// of the authors and should not be interpreted as representing official policies, 
//   either expressed or implied, of the FreeBSD Project.
  
#include <iostream>
#include <fstream>
using namespace std;

#include <Eigen/Core>

#include "epnp.h"
#include "epnp_eigen_debug_tool.h"

const double uc = 599.9;
const double vc = 396.508;
const double fu = 5059.9;
const double fv = 5055.1;


// MtM takes more time than 12x12 opencv SVD with about 180 points and more:

int n = 4;
const double noise = 10;

double rand(double min, double max)
{
  return min + (max - min) * double(rand()) / RAND_MAX;
}

void random_pose(double R[3][3], double t[3])
{
  cout << "random_pose." << endl;
  const double range = 1;

  double phi   = rand(0, range * 3.14159 * 2);
  double theta = rand(0, range * 3.14159);
  double psi   = rand(0, range * 3.14159 * 2);

  R[0][0] = cos(psi) * cos(phi) - cos(theta) * sin(phi) * sin(psi);
  R[0][1] = cos(psi) * sin(phi) + cos(theta) * cos(phi) * sin(psi);
  R[0][2] = sin(psi) * sin(theta);

  R[1][0] = -sin(psi) * cos(phi) - cos(theta) * sin(phi) * cos(psi);
  R[1][1] = -sin(psi) * sin(phi) + cos(theta) * cos(phi) * cos(psi);
  R[1][2] = cos(psi) * sin(theta);

  R[2][0] = sin(theta) * sin(phi);
  R[2][1] = -sin(theta) * cos(phi);
  R[2][2] = cos(theta);

  t[0] = 0.0f;
  t[1] = 0.0f;
  t[2] = 6.0f;
}

void random_point(double & Xw, double & Yw, double & Zw)
{
  double theta = rand(0, 3.14159), phi = rand(0, 2 * 3.14159), R = rand(0, +2);

  Xw =  sin(theta) * sin(phi) * R;
  Yw = -sin(theta) * cos(phi) * R;
  Zw =  cos(theta) * R;
}

void project_with_noise(double R[3][3], double t[3],
			double Xw, double Yw, double Zw,
			double & u, double & v)
{
  double Xc = R[0][0] * Xw + R[0][1] * Yw + R[0][2] * Zw + t[0];
  double Yc = R[1][0] * Xw + R[1][1] * Yw + R[1][2] * Zw + t[1];
  double Zc = R[2][0] * Xw + R[2][1] * Yw + R[2][2] * Zw + t[2];

  double nu = rand(-noise, +noise);
  double nv = rand(-noise, +noise);
  // u = uc + fu * Xc / Zc + nu;
  // v = vc + fv * Yc / Zc + nv;
  u = uc + fu * Xc / Zc ;
  v = vc + fv * Yc / Zc ;
}
// 手动设置3d、2d点
int valuetest(){
  epnp PnP;
  srand(time(0));

  EPnPEigenDebugTool eigen_debug_tool;
  // 读取数值
  Eigen::MatrixXd Pw(n, 3);
  Eigen::MatrixXd pc(n, 2);

  eigen_debug_tool.readFromCSVFile("reference_3d_points.csv", Pw);
  eigen_debug_tool.readFromCSVFile("reference_2d_points.csv", pc);
  // n = Pw.rows();

  PnP.set_internal_parameters(uc, vc, fu, fv);
  PnP.set_maximum_number_of_correspondences(n);

  for(int i = 0; i < n; i++) {
    double Xw, Yw, Zw, u, v;
    Xw = Pw(i,0);
    Yw = Pw(i,1);
    Zw = Pw(i,2);

    u = pc(i,0);
    v = pc(i,1);
    PnP.add_correspondence(Xw, Yw, Zw, u, v);
    cout << " Xw: "<< Xw<< " Yw: "<< Yw << " ZW: " <<Zw <<  " u: " << u << " v: " << v << endl;
  }

  double R_est[3][3], t_est[3];
  double err2 = PnP.compute_pose(R_est, t_est);
  double rot_err, transl_err;

  // PnP.relative_error(rot_err, transl_err, R_true, t_true, R_est, t_est);
  cout << ">>> Reprojection error: " << err2 << endl;
  // cout << ">>> rot_err: " << rot_err << ", transl_err: " << transl_err << endl;
  // cout << endl;
  // cout << "'True reprojection error':"
  //      << PnP.reprojection_error(R_true, t_true) << endl;
  // cout << endl;
  // cout << "True pose:" << endl;
  // PnP.print_pose(R_true, t_true);
  // cout << endl;
  cout << "Found pose:" << endl;
  PnP.print_pose(R_est, t_est);

  
  return 0;

    
}

int randtest(bool setRt,bool setpoint){
  epnp PnP;
  srand(time(0));

  PnP.set_internal_parameters(uc, vc, fu, fv);
  PnP.set_maximum_number_of_correspondences(n);

  double R_true[3][3], t_true[3];
  if(setRt){
    // set R t 
    cout << "****************seting R t :****************"  << endl;
    Eigen::Matrix3d R_true_eigen;
    // R_true_eigen << 0.00143084, 0.0790041, -0.996873, 
    //   0.0860956, 0.993163, 0.0788337 ,
    //   0.996286, -0.0859392, -0.00538085 ;
    cout << "R -R_true_eigen:" << endl;
//     R_true_eigen << 0.0619137, -0.998082, 9.38567e-14, 
// -0.0778879, -0.0048316, -0.99695, 
// 0.995038, 0.0617249, -0.0780376 ;
    R_true_eigen << 0.419437, -0.905384, 0.0659737,
-0.791192, -0.400231, -0.462418,
 0.445071,  0.141757, -0.884204;

    cout << R_true_eigen << endl;
    Eigen::Vector3f t_true_eigen;
    // t_true_eigen << 0.061 ,7.26799, -5.0449;
    // t_true_eigen << -2.06,-13.2623,169.429;
    t_true_eigen << 0, 0, 6;
    
    cout << "t-t_true_eigen:" << endl;
    cout << t_true_eigen << endl;
    for(int i=0; i<R_true_eigen.rows(); i++){
      for(int j=0; j<R_true_eigen.cols(); j++){
        R_true[i][j] = R_true_eigen(i,j);
      }
    }
    for(int i=0; i<t_true_eigen.size(); i++){
      t_true[i] = t_true_eigen(i);
    }
    cout << "R-R_true:" << endl;
    for(int i=0; i<R_true_eigen.rows(); i++){
      for(int j=0; j<R_true_eigen.cols(); j++){
        cout << R_true[i][j] << " ";
      }
      cout << endl;
    }
    cout << "t-t_true:" << endl;
    for(auto j : t_true) cout << j ;
    cout << endl;
  }else
  {
    // use random R t
    random_pose(R_true, t_true);
  }
  // 

  if(setpoint)
  {
    // 读取数值
    EPnPEigenDebugTool eigen_debug_tool;
    Eigen::MatrixXd Pw(4, 3);
    cout << "****************reading point: ****************"<< endl;
    eigen_debug_tool.readFromCSVFile("reference_3d_points.csv", Pw);

    PnP.reset_correspondences();
    cout << "****************seting point: ****************"<< endl;
    n = Pw.rows();
    for(int i = 0; i < n; i++) {
      double Xw, Yw, Zw, u, v;

      // set point or random_point
      Xw = Pw(i,0);
      Yw = Pw(i,1);
      Zw = Pw(i,2);

      project_with_noise(R_true, t_true, Xw, Yw, Zw, u, v);
      cout << " Xw: "<< Xw<< " Yw: "<< Yw << " ZW: " <<Zw <<  " u: " << u << " v: " << v << endl;
      PnP.add_correspondence(Xw, Yw, Zw, u, v);
    }
  }else
  {
    // random point
    cout << "*************random point******" << endl;
    PnP.reset_correspondences();
    for(int i = 0; i < n; i++) {
      double Xw, Yw, Zw, u, v;

      random_point(Xw, Yw, Zw);
      project_with_noise(R_true, t_true, Xw, Yw, Zw, u, v);
      PnP.add_correspondence(Xw, Yw, Zw, u, v);
      cout << "pws: " << Xw << " "<<  Yw << "  "<< Zw<< endl;
      cout << "us: " << u <<" "<< v<< endl;
    }
  }
  

  double R_est[3][3], t_est[3];
  double err2 = PnP.compute_pose(R_est, t_est);
  double rot_err, transl_err;

  PnP.relative_error(rot_err, transl_err, R_true, t_true, R_est, t_est);
  cout << ">>> Reprojection error: " << err2 << endl;
  cout << ">>> rot_err: " << rot_err << ", transl_err: " << transl_err << endl;
  cout << endl;
  cout << "'True reprojection error':"
       << PnP.reprojection_error(R_true, t_true) << endl;
  cout << endl;
  cout << "True pose:" << endl;
  PnP.print_pose(R_true, t_true);
  cout << endl;
  cout << "Found pose:" << endl;
  PnP.print_pose(R_est, t_est);

  return 0;
}
int main(int argc, char ** argv, char** isrand  ,char**  issetRt , char**  issetpoint )
{
  bool bvaluetest,bsetRt, bsetpoint;
  cout << "argv[1]" << argv[1] << argc<< endl;
  if(argc==2 && string(argv[1]) == "true")
  {
    bvaluetest  = true;
  }
  else if(argc!=4)
  {
    bvaluetest  = false ;
    bsetRt = false;
    bsetpoint =  false;
  }else{
    std::stringstream srand(argv[1]);
    std::stringstream sRt(argv[2]);
    std::stringstream spoint(argv[3]);
    srand >> std::boolalpha >> bvaluetest ;
    sRt >> std::boolalpha >> bsetRt ;
    spoint >>  std::boolalpha >>  bsetpoint;
  }  
        

  std::cout << "--args: " << std::endl   
      << "  bvaluetest , 是否使用已有数据测试(true:value|false:random): " << bvaluetest << std::endl 
      << "  setRt , 是否设置变换矩阵： " << bsetRt << std::endl
      << "  setPoint, 是否设置3d点坐标: " << bsetpoint << std::endl;

  if(bvaluetest){
    return valuetest();
  }else{
    return randtest(bsetRt,bsetpoint);
  }
}

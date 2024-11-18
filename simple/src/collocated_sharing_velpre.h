#ifndef COLLOCATED_SHARING_VELPRE_H
#define COLLOCATED_SHARING_VELPRE_H

#include <vector>

#include "StructuredMesh.h"


/**
 * 设置面速度的边界条件
 * @return c                 根据边界条件，对面速度进行设置
 */
void set_face_vel_bc(StructuredMesh& mesh, FluidBoundary& fluidBoundary, int itest, int dim, int ncx, int ncy, int ncz,
                     const std::vector<std::vector<std::vector<double>>>& u,
                     const std::vector<std::vector<std::vector<double>>>& v,
                     const std::vector<std::vector<std::vector<double>>>& w,
                     std::vector<std::vector<std::vector<double>>>& uf,
                     std::vector<std::vector<std::vector<double>>>& vf,
                     std::vector<std::vector<std::vector<double>>>& wf);


/**
 * 计算动量方程的系数
 * @return c                 离散动量方程并计算系数矩阵，用于数值求解速度场
 */
void momentum_coefs(StructuredMesh& caseData, int dim, int ncx, int ncy, int ncz, int ncoef, double dt,
                    const std::vector<std::vector<std::vector<double>>>& mu,
                    const std::vector<std::vector<std::vector<double>>>& dens,
                    const std::vector<std::vector<std::vector<double>>>& u,
                    const std::vector<std::vector<std::vector<double>>>& v,
                    const std::vector<std::vector<std::vector<double>>>& w,
                    const std::vector<std::vector<std::vector<double>>>& uf,
                    const std::vector<std::vector<std::vector<double>>>& vf,
                    const std::vector<std::vector<std::vector<double>>>& wf,
                    std::vector<std::vector<std::vector<std::vector<double>>>>& cu,
                    std::vector<std::vector<std::vector<std::vector<double>>>>& cv,
                    std::vector<std::vector<std::vector<std::vector<double>>>>& cw);



/**
 * 根据边界条件修正动量方程的系数矩阵
 */ 
void momentum_coef_bcs(StructuredMesh& caseData, FluidBoundary& fluidBoundary, int itest, int dim, int ncx, int ncy, int ncz, int dir,
                       std::vector<std::vector<std::vector<std::vector<double>>>>& cu);


/**
 * 计算动量方程中因压力梯度引起的源项
 */
void momentum_gradp(StructuredMesh& caseData, FluidBoundary& fluidBoundary, int itest, int dim, int ncx, int ncy, int ncz,
                    const std::vector<std::vector<std::vector<double>>>& p,
                    std::vector<std::vector<std::vector<std::vector<double>>>>& cu,
                    std::vector<std::vector<std::vector<std::vector<double>>>>& cv,
                    std::vector<std::vector<std::vector<std::vector<double>>>>& cw);


/**
 * 计算Rhie-Chow插值的面速度
 */
void rhie_chow_face_velocity(StructuredMesh& caseData, FluidBoundary& fluidBoundary, int itest, int dim, int ncx, int ncy, int ncz, int ncoef,
                             const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, double dt,
                             const std::vector<std::vector<std::vector<double>>>& dens,
                             const std::vector<std::vector<std::vector<double>>>& p,
                             const std::vector<std::vector<std::vector<double>>>& u, const std::vector<std::vector<std::vector<double>>>& v,
                             const std::vector<std::vector<std::vector<double>>>& w,
                             const std::vector<std::vector<std::vector<double>>>& u0, const std::vector<std::vector<std::vector<double>>>& v0,
                             const std::vector<std::vector<std::vector<double>>>& w0,
                             std::vector<std::vector<std::vector<double>>>& uf, std::vector<std::vector<std::vector<double>>>& vf,
                             std::vector<std::vector<std::vector<double>>>& wf,
                             const std::vector<std::vector<std::vector<std::vector<double>>>>& cu,
                             const std::vector<std::vector<std::vector<std::vector<double>>>>& cv,
                             const std::vector<std::vector<std::vector<std::vector<double>>>>& cw);

/**
 * 离散压力校正方程并计算其系数，用于压力场修正
 */
void pressure_coefs(StructuredMesh& caseData, FluidBoundary& fluidBoundary, int itest, int dim, int ncx, int ncy, int ncz,
                    const std::vector<std::vector<std::vector<double>>>& dens,
                    const std::vector<std::vector<std::vector<double>>>& uf,
                    const std::vector<std::vector<std::vector<double>>>& vf,
                    const std::vector<std::vector<std::vector<double>>>& wf,
                    const std::vector<std::vector<std::vector<std::vector<double>>>>& cu,
                    const std::vector<std::vector<std::vector<std::vector<double>>>>& cv,
                    const std::vector<std::vector<std::vector<std::vector<double>>>>& cw,
                    std::vector<std::vector<std::vector<std::vector<double>>>>& cp);


/**
 * 使用压力校正结果修正压力场
 */
void correct_pressure(int ncx, int ncy, int ncz, double relax_p,
                      std::vector<std::vector<std::vector<double>>>& pp,
                      std::vector<std::vector<std::vector<double>>>& p);

/**
 * 使用压力校正结果修正速度场
 */
void correct_velocity(int dim, int ncx, int ncy, int ncz,
                      const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z,
                      const std::vector<std::vector<std::vector<std::vector<double>>>>& cu,
                      const std::vector<std::vector<std::vector<std::vector<double>>>>& cv,
                      const std::vector<std::vector<std::vector<double>>>& pp,
                      std::vector<std::vector<std::vector<double>>>& u,
                      std::vector<std::vector<std::vector<double>>>& v,
                      std::vector<std::vector<std::vector<double>>>& uf,
                      std::vector<std::vector<std::vector<double>>>& vf);


#endif 

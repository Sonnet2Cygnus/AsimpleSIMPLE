#ifndef SOLVERS_H
#define SOLVERS_H

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>

#include "StructuredMesh.h"


/**
 * 计算标量场的 L2 范数，并更新 case 对象中的 l2_u 和 l2_max_u
 * @param caseData          网格数据对象
 * @param it_nl             非线性迭代次数
 * @param ncx, ncy, ncz     网格单元数
 * @param u0                前一步的标量场
 * @param u                 当前的标量场
 * @return                  返回一个包含当前 L2 范数和最大 L2 范数的 pair
 */
std::pair<double, double> eqn_scalar_norm2(StructuredMesh& caseData, int dim, int it_nl, int ncx, int ncy, int ncz,
                                           const std::vector<std::vector<std::vector<double>>>& u0,
                                           const std::vector<std::vector<std::vector<double>>>& u,
                                           const std::string& var);



/**
 * 计算压力项的 L2 范数
 * @param c         系数矩阵
 */
void eqn_scalar_norm2_cp(StructuredMesh& caseData, int dim, int it_nl, int ncx, int ncy, int ncz,
                         const std::vector<std::vector<std::vector<std::vector<double>>>>& c,
                         const std::string& var);

/**
 * 使用雅可比方法求解标量方程
 * @param niter     最大迭代次数
 * @param relax     松弛因子
 * @param ncoef     系数矩阵的列数
 * @param ct        系数矩阵
 * @param t         标量场
 * @param initzero  是否初始化为零
 * @param res       容差
 */
void scalar_pj(StructuredMesh& caseData, int dim, int it_nl, int niter, double relax,
               int ncx, int ncy, int ncz, int ncoef,
               const std::vector<std::vector<std::vector<std::vector<double>>>>& ct,
               std::vector<std::vector<std::vector<double>>>& t,
               bool initzero, double res);


// 使用高斯-赛德尔方法求解标量方程
void scalar_gs(StructuredMesh& caseData, int dim, int it_nl, int niter, double relax,
               int ncx, int ncy, int ncz, int ncoef,
               const std::vector<std::vector<std::vector<std::vector<double>>>>& ct,
               std::vector<std::vector<std::vector<double>>>& t,
               bool initzero, double res);

#endif

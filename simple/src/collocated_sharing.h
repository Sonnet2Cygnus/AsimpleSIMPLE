#ifndef COLLOCATED_SHARING_H
#define COLLOCATED_SHARING_H

#include <vector>

#include "StructuredMesh.h"
#include "kernels.h"

/**
 * 计算传导系数
 * @param spht               比热容
 * @param con                导热系数
 * @param heat_src           热源项
 * @param dens                密度场
 * @return ct                  计算出的传导系数矩阵，包含每个网格单元的传导系数
 */
void conduction_coefs(StructuredMesh& caseData, int dim, int ncx, int ncy, int ncz, int ncoef, double dt,
                      double spht, double con, double heat_src,
                      const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z,
                      const std::vector<std::vector<std::vector<double>>>& dens,
                      const std::vector<std::vector<std::vector<double>>>& t,
                      const std::vector<std::vector<std::vector<double>>>& uf,
                      const std::vector<std::vector<std::vector<double>>>& vf,
                      const std::vector<std::vector<std::vector<double>>>& wf,
                      std::vector<std::vector<std::vector<std::vector<double>>>>& ct);

/**
 * 对传导系数矩阵应用边界条件
 * @return ct                  应用边界条件后的传导系数矩阵
 */
void conduction_coef_bcs(StructuredMesh& caseData, FluidBoundary& fluidboundary, int dim, int ncx, int ncy, int ncz, int ncoef,
                         double dt, double con,
                         const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z,
                         const std::vector<std::vector<std::vector<double>>>& dens,
                         const std::vector<std::vector<std::vector<double>>>& t,
                         std::vector<std::vector<std::vector<std::vector<double>>>>& ct);

#endif
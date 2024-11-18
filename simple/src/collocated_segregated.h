#ifndef COLLOCATED_SEGREGATED_H
#define COLLOCATED_SEGREGATED_H

#include "StructuredMesh.h"
#include "collocated_sharing.h"
#include "collocated_sharing_velpre.h"
#include "Solvers.h"
#include "kernels.h"


/**
 * @brief 主控函数，用于执行基于共置网格的分离算法
 *
 * @param it                    当前的迭代次数
 * @param caseData 当           前计算用例的网格数据
 * @param fluid                 流体属性，包括密度、粘性等参数
 * @param fluidBoundary         流体边界条件，包括速度、压力、温度等边界数据
 *
 * @return bool 是否成功完成当前步骤的计算
 *
 * 该函数为基于共置网格的分离算法的主流程，包含以下部分：
 * 1. 动量方程的离散与求解
 * 2. 压力校正方程的离散与求解
 * 3. 速度场和压力场的校正
 * 4. 满足连续性条件的校验
 */
bool collocated_segregated(int it, StructuredMesh& caseData, Fluid& fluid, FluidBoundary& fluidBoundary);

#endif 

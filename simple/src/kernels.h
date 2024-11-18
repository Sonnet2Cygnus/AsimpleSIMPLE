#ifndef KERNELS_H
#define KERNELS_H

#include <algorithm>
#include <cmath>    

/**
 * 根据网格尺寸 dx, dy, dz 计算网格面的面积和体积
 * @return                  返回 x,y,z面积和网格体积
 */
void cal_area_vol(double dx, double dy, double dz,
                  double& area_x, double& area_y, double& area_z, double& vol);

/**
 * 根据 Peclet 数计算修正系数 ap
 * @param pec                Peclet 数，用于描述对流和扩散的相对强度
 * @return                  修正系数 ap
 */
double a_pec_pow(double pec);


/**
 * 计算相邻网格单元的系数 a
 * @param idx                面积法向上的倒数 如 1/dx
 * @param ul                 左侧网格的速度
 * @param gl                 左侧网格的扩散系数
 * @param sign_f             流动方向的符号
 * @return                  邻单元的修正系数 a
 */
double a_nb(int conv_scheme, double area, double idx,
            double ul, double ur, double gl, double gr, double rho, double sign_f);

#endif 

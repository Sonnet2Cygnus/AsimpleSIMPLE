#include "kernels.h"

void cal_area_vol(double dx, double dy, double dz,
                  double& area_x, double& area_y, double& area_z, double& vol) {
    area_x = dy * dz;
    area_y = dx * dz;
    area_z = dx * dy;
    vol    = dx * area_x;
}


double a_pec_pow(double pec) {

    // 线性削弱系数
    double ap = 1.0 - 0.1 * std::abs(pec);

    // 保证非负 5次幂
    ap = std::max(0.0, std::pow(ap, 5));

    return ap;
}

double a_nb(int conv_scheme, double area, double idx,
            double ul, double ur, double gl, double gr, double rho, double sign_f) {

    // 计算质量流量
    double f = rho * 0.5 * (ul + ur);

    // 扩散系数的谐调平均，避免零除
    double d = 2.0 * gl * gr / (gl + gr + 1.0e-12) * idx;

    // 修正系数
    double a = 0.0;

    // 根据对流格式选择不同的计算方式 0-迎风 1-线性 2-修正 3- 占位
    if (conv_scheme == 0) {
        a = area * (d + std::max(0.0, sign_f * f));
    } else if (conv_scheme == 1) {
        a = area * (d * (1.0 - 0.5 * std::abs(f / d)) + std::max(0.0, sign_f * f));
    } else if (conv_scheme == 2) {
        a = area * (d * a_pec_pow(std::abs(f / d)) + std::max(0.0, sign_f * f));
    } else if (conv_scheme == 3) {
        a = 0.0; 
    }

    return a;
}

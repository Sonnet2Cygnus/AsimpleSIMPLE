#include "collocated_sharing.h"
#include <cmath>
#include <algorithm>

void conduction_coefs(StructuredMesh& caseData, int dim, int ncx, int ncy, int ncz, int ncoef, double dt,
                      double spht, double con, double heat_src,
                      const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z,
                      const std::vector<std::vector<std::vector<double>>>& dens,
                      const std::vector<std::vector<std::vector<double>>>& t,
                      const std::vector<std::vector<std::vector<double>>>& uf,
                      const std::vector<std::vector<std::vector<double>>>& vf,
                      const std::vector<std::vector<std::vector<double>>>& wf,
                      std::vector<std::vector<std::vector<std::vector<double>>>>& ct) {

    int id_aP = caseData.id_aP;
    int id_aE = caseData.id_aE;
    int id_aW = caseData.id_aW;
    int id_aN = caseData.id_aN;
    int id_aS = caseData.id_aS;
    int id_aT = 0;
    int id_aB = 0;
    if (dim == 3) {
        id_aT = caseData.id_aT;
        id_aB = caseData.id_aB;
    }
    int id_bsrc = caseData.id_bsrc;

    int conv_scheme = caseData.conv_scheme;  

    // 时间步长的倒数
    double idt = 1.0 / dt;

    // 计算网格尺寸和几何信息
    double dx = x[1] - x[0];
    double dy = y[1] - y[0];
    double dz = z[1] - z[0];

    double area_x, area_y, area_z, vol;
    cal_area_vol(dx, dy, dz, area_x, area_y, area_z, vol);

    double idx = 1.0 / dx;
    double idy = 1.0 / dy;
    double idz = 1.0 / dz;

    // 遍历所有网格单元
    for (int k = 0; k < ncz; ++k) {
        for (int j = 0; j < ncy; ++j) {
            for (int i = 0; i < ncx; ++i) {

                double aE = 0.0, aW = 0.0;
                double aN = 0.0, aS = 0.0;
                double aT = 0.0, aB = 0.0;
                double aP = 0.0;

                double sC = heat_src;       // 热源项
                double sP = 0.0;            // 线性化源项
                double bsrc = 0.0;

                double rho = 0.0;
                double mul = con, mur = con;
                double ul = 0.0, ur = 0.0;

                // 东侧面
                if (i == ncx - 1) {
                    rho = dens[i][j][k];
                } else {
                    rho = 0.5 * (dens[i][j][k] + dens[i + 1][j][k]);
                }
                ul = uf[i + 1][j][k];
                ur = uf[i + 1][j][k];
                aE = a_nb(conv_scheme, area_x, idx, ul, ur, mul, mur, rho, -1.0);

                // 西
                if (i == 0) {
                    rho = dens[i][j][k];
                } else {
                    rho = 0.5 * (dens[i][j][k] + dens[i - 1][j][k]);
                }
                ul = uf[i][j][k];
                ur = uf[i][j][k];
                aW = a_nb(conv_scheme, area_x, idx, ul, ur, mul, mur, rho, 1.0);

                // 北
                if (j == ncy - 1) {
                    rho = dens[i][j][k];
                } else {
                    rho = 0.5 * (dens[i][j][k] + dens[i][j + 1][k]);
                }
                ul = vf[i][j + 1][k];
                ur = vf[i][j + 1][k];
                aN = a_nb(conv_scheme, area_y, idy, ul, ur, mul, mur, rho, -1.0);

                // 南
                if (j == 0) {
                    rho = dens[i][j][k];
                } else {
                    rho = 0.5 * (dens[i][j][k] + dens[i][j - 1][k]);
                }
                ul = vf[i][j][k];
                ur = vf[i][j][k];
                aS = a_nb(conv_scheme, area_y, idy, ul, ur, mul, mur, rho, 1.0);

                if (dim == 3) {
                    // 顶
                    if (k == ncz - 1) {
                        rho = dens[i][j][k];
                    } else {
                        rho = 0.5 * (dens[i][j][k] + dens[i][j][k + 1]);
                    }
                    ul = wf[i][j][k + 1];
                    ur = wf[i][j][k + 1];
                    aT = a_nb(conv_scheme, area_z, idz, ul, ur, mul, mur, rho, -1.0);

                    // 底
                    if (k == 0) {
                        rho = dens[i][j][k];
                    } else {
                        rho = 0.5 * (dens[i][j][k] + dens[i][j][k - 1]);
                    }
                    ul = wf[i][j][k];
                    ur = wf[i][j][k];
                    aB = a_nb(conv_scheme, area_z, idz, ul, ur, mul, mur, rho, 1.0);
                }

                // 计算主对角线系数
                rho = dens[i][j][k];
                double aP0 = rho * spht * vol * idt;                    // 时间导数项
                aP = aE + aW + aN + aS + aT + aB + aP0 - sP * vol;

                // 存储系数到系数矩阵
                ct[i][j][k][id_aP] = aP;
                ct[i][j][k][id_aE] = -aE;
                ct[i][j][k][id_aW] = -aW;
                ct[i][j][k][id_aN] = -aN;
                ct[i][j][k][id_aS] = -aS;
                if (dim == 3) {
                    ct[i][j][k][id_aT] = -aT;
                    ct[i][j][k][id_aB] = -aB;
                }

                // 计算源项
                bsrc = sC * vol + aP0 * t[i][j][k];
                ct[i][j][k][id_bsrc] = bsrc;
            }
        }
    }
}

void conduction_coef_bcs(StructuredMesh& caseData, FluidBoundary& fluidboundary, int dim, int ncx, int ncy, int ncz, int ncoef,
                         double dt, double con,
                         const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z,
                         const std::vector<std::vector<std::vector<double>>>& dens,
                         const std::vector<std::vector<std::vector<double>>>& t,
                         std::vector<std::vector<std::vector<std::vector<double>>>>& ct) {

    int id_aP = caseData.id_aP;
    int id_aE = caseData.id_aE;
    int id_aW = caseData.id_aW;
    int id_aN = caseData.id_aN;
    int id_aS = caseData.id_aS;
    int id_aT = 0;
    int id_aB = 0;
    
    if (dim == 3) {
        id_aT = caseData.id_aT;
        id_aB = caseData.id_aB;
    }

    int id_bsrc = caseData.id_bsrc;     // 源项

    // 网格中心坐标
    const std::vector<double>& xc = caseData.xc;
    const std::vector<double>& yc = caseData.yc;
    const std::vector<double>& zc = caseData.zc;

    // 边界标识以及各个边界
    const auto& bcid = fluidboundary.bcid;
    int fid_e = fluidboundary.fid_e;
    int fid_w = fluidboundary.fid_w;
    int fid_n = fluidboundary.fid_n;
    int fid_s = fluidboundary.fid_s;
    int fid_t = fluidboundary.fid_t;
    int fid_b = fluidboundary.fid_b;

    // 边界条件类型
    const auto& bcs = fluidboundary.bcs;
    const auto& bcs_temp = fluidboundary.bcs_temp;
    int bc_none = fluidboundary.bc_none;
    int bc_wall = fluidboundary.bc_wall;
    int bc_inlet = fluidboundary.bc_inlet;
    int bc_outlet = fluidboundary.bc_outlet;
    int temp_bc_constant = fluidboundary.bc_constant;
    int temp_bc_heatflux = fluidboundary.bc_flux;

    // 遍历所有网格单元
    for (int k = 0; k < ncz; ++k) {
        for (int j = 0; j < ncy; ++j) {
            for (int i = 0; i < ncx; ++i) {

                // 获取网格单元的边界条件标识
                int bcid_e = bcid[i][j][k][fid_e];
                int bcid_w = bcid[i][j][k][fid_w];
                int bcid_n = bcid[i][j][k][fid_n];
                int bcid_s = bcid[i][j][k][fid_s];

                // 获取边界类型
                int bc_e = bcs[bcid_e].type;
                int bc_w = bcs[bcid_w].type;
                int bc_n = bcs[bcid_n].type;
                int bc_s = bcs[bcid_s].type;

                int ttype_e = bcs_temp[bcid_e].temp_type;
                int ttype_w = bcs_temp[bcid_w].temp_type;
                int ttype_n = bcs_temp[bcid_n].temp_type;
                int ttype_s = bcs_temp[bcid_s].temp_type;

                // 获取边界值
                double bc_te = bcs_temp[bcid_e].t;
                double bc_tw = bcs_temp[bcid_w].t;
                double bc_tn = bcs_temp[bcid_n].t;
                double bc_ts = bcs_temp[bcid_s].t;

                double bc_fluxe = bcs_temp[bcid_e].heat_flux;
                double bc_fluxw = bcs_temp[bcid_w].heat_flux;
                double bc_fluxn = bcs_temp[bcid_n].heat_flux;
                double bc_fluxs = bcs_temp[bcid_s].heat_flux;

                // 东侧边界条件
                if (bc_e == bc_inlet) {
                    ct[i][j][k][id_aP] -= ct[i][j][k][id_aE];
                    ct[i][j][k][id_bsrc] -= 2.0 * ct[i][j][k][id_aE] * bc_te;
                    ct[i][j][k][id_aE] = 0.0;
                } else if (bc_e == bc_outlet) {
                    ct[i][j][k][id_aP] += ct[i][j][k][id_aE];
                    ct[i][j][k][id_aE] = 0.0;
                } else if (bc_e == bc_wall) {
                    if (ttype_e == temp_bc_constant) {
                        ct[i][j][k][id_aP] -= ct[i][j][k][id_aE];
                        ct[i][j][k][id_bsrc] -= 2.0 * ct[i][j][k][id_aE] * bc_te;
                        ct[i][j][k][id_aE] = 0.0;
                    } else if (ttype_e == temp_bc_heatflux) {
                        ct[i][j][k][id_aP] += ct[i][j][k][id_aE];
                        ct[i][j][k][id_bsrc] += ct[i][j][k][id_aE] * bc_fluxe * 2.0 * (x[i + 1] - xc[i]) / con;
                        ct[i][j][k][id_aE] = 0.0;
                    }
                }

                // 西
                if (bc_w == bc_inlet) {
                    ct[i][j][k][id_aP] -= ct[i][j][k][id_aW];
                    ct[i][j][k][id_bsrc] -= 2.0 * ct[i][j][k][id_aW] * bc_tw;
                    ct[i][j][k][id_aW] = 0.0;
                } else if (bc_w == bc_outlet) {
                    ct[i][j][k][id_aP] += ct[i][j][k][id_aW];
                    ct[i][j][k][id_aW] = 0.0;
                } else if (bc_w == bc_wall) {
                    if (ttype_w == temp_bc_constant) {
                        ct[i][j][k][id_aP] -= ct[i][j][k][id_aW];
                        ct[i][j][k][id_bsrc] -= 2.0 * ct[i][j][k][id_aW] * bc_tw;
                        ct[i][j][k][id_aW] = 0.0;
                    } else if (ttype_w == temp_bc_heatflux) {
                        ct[i][j][k][id_aP] += ct[i][j][k][id_aW];
                        ct[i][j][k][id_bsrc] += ct[i][j][k][id_aW] * bc_fluxw * 2.0 * (x[i] - xc[i]) / con;
                        ct[i][j][k][id_aW] = 0.0;
                    }
                }

                // 北
                if (bc_n == bc_inlet) {
                    ct[i][j][k][id_aP] -= ct[i][j][k][id_aN];
                    ct[i][j][k][id_bsrc] -= 2.0 * ct[i][j][k][id_aN] * bc_tn;
                    ct[i][j][k][id_aN] = 0.0;
                } else if (bc_n == bc_outlet) {
                    ct[i][j][k][id_aP] += ct[i][j][k][id_aN];
                    ct[i][j][k][id_aN] = 0.0;
                } else if (bc_n == bc_wall) {
                    if (ttype_n == temp_bc_constant) {
                        ct[i][j][k][id_aP] -= ct[i][j][k][id_aN];
                        ct[i][j][k][id_bsrc] -= 2.0 * ct[i][j][k][id_aN] * bc_tn;
                        ct[i][j][k][id_aN] = 0.0;
                    } else if (ttype_n == temp_bc_heatflux) {
                        ct[i][j][k][id_aP] += ct[i][j][k][id_aN];
                        ct[i][j][k][id_bsrc] += ct[i][j][k][id_aN] * bc_fluxn * 2.0 * (y[j + 1] - yc[j]) / con;
                        ct[i][j][k][id_aN] = 0.0;
                    }
                }

                // 南
                if (bc_s == bc_inlet) {
                    ct[i][j][k][id_aP] -= ct[i][j][k][id_aS];
                    ct[i][j][k][id_bsrc] -= 2.0 * ct[i][j][k][id_aS] * bc_ts;
                    ct[i][j][k][id_aS] = 0.0;
                } else if (bc_s == bc_outlet) {
                    ct[i][j][k][id_aP] += ct[i][j][k][id_aS];
                    ct[i][j][k][id_aS] = 0.0;
                } else if (bc_s == bc_wall) {
                    if (ttype_s == temp_bc_constant) {
                        ct[i][j][k][id_aP] -= ct[i][j][k][id_aS];
                        ct[i][j][k][id_bsrc] -= 2.0 * ct[i][j][k][id_aS] * bc_ts;
                        ct[i][j][k][id_aS] = 0.0;
                    } else if (ttype_s == temp_bc_heatflux) {
                        ct[i][j][k][id_aP] += ct[i][j][k][id_aS];
                        ct[i][j][k][id_bsrc] += ct[i][j][k][id_aS] * bc_fluxs * 2.0 * (y[j] - yc[j]) / con;
                        ct[i][j][k][id_aS] = 0.0;
                    }
                }
                // 3D 的暂时保留
            }
        }
    }
}

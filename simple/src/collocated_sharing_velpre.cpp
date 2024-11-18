#include <cmath>

#include "collocated_sharing_velpre.h"


void set_face_vel_bc(StructuredMesh& caseData, FluidBoundary& fluidBoundary, int itest, int dim, int ncx, int ncy, int ncz,
                     const std::vector<std::vector<std::vector<double>>>& u,
                     const std::vector<std::vector<std::vector<double>>>& v,
                     const std::vector<std::vector<std::vector<double>>>& w,
                     std::vector<std::vector<std::vector<double>>>& uf,
                     std::vector<std::vector<std::vector<double>>>& vf,
                     std::vector<std::vector<std::vector<double>>>& wf) {

    const auto& bcid = fluidBoundary.bcid;

    int fid_e = fluidBoundary.fid_e; // 冬/右
    int fid_w = fluidBoundary.fid_w; 
    int fid_n = fluidBoundary.fid_n; 
    int fid_s = fluidBoundary.fid_s; 
    int fid_t = 0;
    int fid_b = 0;
    if (dim == 3) {
        fid_t = fluidBoundary.fid_t; 
        fid_b = fluidBoundary.fid_b; 
    }

    const auto& bcs = fluidBoundary.bcs;                    // 边界条件
    const auto& bcs_vel = fluidBoundary.bcs_vel;            // 速度边界条件

    int bc_none = fluidBoundary.bc_none;   
    int bc_wall = fluidBoundary.bc_wall;   
    int bc_inlet = fluidBoundary.bc_inlet; 
    int bc_outlet = fluidBoundary.bc_outlet; 

    // 质量流量
    double mf_in = 0.0;
    double mf_out = 0.0;

    // 遍历
    for (int k = 0; k < ncz; ++k) {
        for (int j = 0; j < ncy; ++j) {
            for (int i = 0; i < ncx; ++i) {

                // 编号
                int bcid_e = bcid[i][j][k][fid_e]; 
                int bcid_w = bcid[i][j][k][fid_w]; 
                int bcid_n = bcid[i][j][k][fid_n]; 
                int bcid_s = bcid[i][j][k][fid_s]; 

                // 边界
                int bc_e = bcs[bcid_e].type; 
                int bc_w = bcs[bcid_w].type; 
                int bc_n = bcs[bcid_n].type; 
                int bc_s = bcs[bcid_s].type; 

                // 边界速度
                double bc_ue = bcs_vel[bcid_e].u; 
                double bc_uw = bcs_vel[bcid_w].u; 
                double bc_un = bcs_vel[bcid_n].u; 
                double bc_us = bcs_vel[bcid_s].u; 

                // 东侧边界条件
                if (bc_e == bc_inlet) {

                    // 设置入口速度
                    uf[i + 1][j][k] = bc_ue;

                    // 累加进口mf
                    mf_in += uf[i + 1][j][k];
                } else if (bc_e == bc_outlet) {
                    uf[i + 1][j][k] = u[i][j][k];
                    mf_out += uf[i + 1][j][k];
                } else if (bc_e == bc_wall) {
                    
                    // 设定为常数
                    uf[i + 1][j][k] = bc_ue;
                }

                // 西
                if (bc_w == bc_inlet) {
                    uf[i][j][k] = bc_uw;
                    mf_in += uf[i][j][k];
                } else if (bc_w == bc_outlet) {
                    uf[i][j][k] = u[i][j][k];
                    mf_out += uf[i][j][k];
                } else if (bc_w == bc_wall) {
                    uf[i][j][k] = bc_uw;
                }

                // 北
                if (bc_n == bc_inlet) {
                    vf[i][j + 1][k] = bc_un;
                    mf_in += vf[i][j + 1][k];
                } else if (bc_n == bc_outlet) {
                    vf[i][j + 1][k] = v[i][j][k];
                    mf_out += vf[i][j + 1][k];
                } else if (bc_n == bc_wall) {
                    vf[i][j + 1][k] = bc_un;
                }

                // 南
                if (bc_s == bc_inlet) {
                    vf[i][j][k] = bc_us;
                    mf_in += vf[i][j][k];
                } else if (bc_s == bc_outlet) {
                    vf[i][j][k] = v[i][j][k];
                    mf_out += vf[i][j][k];
                } else if (bc_s == bc_wall) {
                    vf[i][j][k] = bc_us;
                }

                
                if (dim == 3) {
                    // 此处留3D
                }
            }
        }
    }

    // 累积进口和出口的质量通量
    for (int k = 0; k < ncz; ++k) {
        for (int j = 0; j < ncy; ++j) {
            mf_in += uf[0][j][k];
            mf_out += uf[ncx][j][k];
        }
    }
}


// 计算邻接单元的系数 a_nb
double a_nb(int conv_scheme, double area, double idh, double ul, double ur, double mul, double mur, double rho, double sign);


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
                    std::vector<std::vector<std::vector<std::vector<double>>>>& cw) {

    int id_aP = caseData.id_aP;
    int id_aE = caseData.id_aE;
    int id_aW = caseData.id_aW;
    int id_aN = caseData.id_aN;
    int id_aS = caseData.id_aS;
    int id_aT = caseData.id_aT;
    int id_aB = caseData.id_aB;
    int id_bsrc = caseData.id_bsrc;

    int conv_scheme = caseData.conv_scheme;

    double idt = 1.0 / dt;

    double dx = caseData.dx;
    double dy = caseData.dy;
    double dz = caseData.dz;

    double area_x = dy * dz;
    double area_y = dx * dz;
    double area_z = dx * dy;
    double vol = dx * dy * dz;

    double idx = 1.0 / dx;
    double idy = 1.0 / dy;
    double idz = 1.0 / dz;

    for (int k = 0; k < ncz; ++k) {
        for (int j = 0; j < ncy; ++j) {
            for (int i = 0; i < ncx; ++i) {

                double aE = 0.0, aW = 0.0;
                double aN = 0.0, aS = 0.0;
                double aT = 0.0, aB = 0.0;
                double aP = 0.0;

                double sC = 0.0;
                double sP = 0.0;
                double bsrc = 0.0;

                double rho, mul, mur;
                double ul, ur;

                // E
                if (i == ncx - 1) {
                    rho = dens[i][j][k];
                    mul = mu[i][j][k];
                    mur = mu[i][j][k];
                } else {
                    rho = 0.5 * (dens[i][j][k] + dens[i + 1][j][k]);
                    mul = mu[i][j][k];
                    mur = mu[i + 1][j][k];
                }

                ul = uf[i + 1][j][k];
                ur = uf[i + 1][j][k];
                aE = a_nb(conv_scheme, area_x, idx, ul, ur, mul, mur, rho, -1.0);

                // W
                if (i == 0) {
                    rho = dens[i][j][k];
                    mul = mu[i][j][k];
                    mur = mu[i][j][k];
                } else {
                    rho = 0.5 * (dens[i][j][k] + dens[i - 1][j][k]);
                    mul = mu[i - 1][j][k];
                    mur = mu[i][j][k];
                }

                ul = uf[i][j][k];
                ur = uf[i][j][k];
                aW = a_nb(conv_scheme, area_x, idx, ul, ur, mul, mur, rho, 1.0);

                // N
                if (j == ncy - 1) {
                    rho = dens[i][j][k];
                    mul = mu[i][j][k];
                    mur = mu[i][j][k];
                } else {
                    rho = 0.5 * (dens[i][j][k] + dens[i][j + 1][k]);
                    mul = mu[i][j][k];
                    mur = mu[i][j + 1][k];
                }

                ul = vf[i][j + 1][k];
                ur = vf[i][j + 1][k];
                aN = a_nb(conv_scheme, area_y, idy, ul, ur, mul, mur, rho, -1.0);

                // S
                if (j == 0) {
                    rho = dens[i][j][k];
                    mul = mu[i][j][k];
                    mur = mu[i][j][k];
                } else {
                    rho = 0.5 * (dens[i][j][k] + dens[i][j - 1][k]);
                    mul = mu[i][j - 1][k];
                    mur = mu[i][j][k];
                }

                ul = vf[i][j][k];
                ur = vf[i][j][k];
                aS = a_nb(conv_scheme, area_y, idy, ul, ur, mul, mur, rho, 1.0);

                if (dim == 3) {
                    // Top face
                    if (k == ncz - 1) {
                        rho = dens[i][j][k];
                        mul = mu[i][j][k];
                        mur = mu[i][j][k];
                    } else {
                        rho = 0.5 * (dens[i][j][k] + dens[i][j][k + 1]);
                        mul = mu[i][j][k];
                        mur = mu[i][j][k + 1];
                    }

                    ul = wf[i][j][k + 1];
                    ur = wf[i][j][k + 1];
                    aT = a_nb(conv_scheme, area_z, idz, ul, ur, mul, mur, rho, -1.0);

                    // Bottom face
                    if (k == 0) {
                        rho = dens[i][j][k];
                        mul = mu[i][j][k];
                        mur = mu[i][j][k];
                    } else {
                        rho = 0.5 * (dens[i][j][k] + dens[i][j][k - 1]);
                        mul = mu[i][j][k - 1];
                        mur = mu[i][j][k];
                    }

                    ul = wf[i][j][k];
                    ur = wf[i][j][k];
                    aB = a_nb(conv_scheme, area_z, idz, ul, ur, mul, mur, rho, 1.0);
                }

                rho = dens[i][j][k];
                double aP0 = rho * vol * idt;
                aP = aE + aW + aN + aS + aT + aB + aP0 - sP * vol;

                // 更新cu矩阵
                cu[i][j][k][id_aP] = aP;
                cu[i][j][k][id_aE] = -aE;
                cu[i][j][k][id_aW] = -aW;
                cu[i][j][k][id_aN] = -aN;
                cu[i][j][k][id_aS] = -aS;
                if (dim == 3) {
                    cu[i][j][k][id_aT] = -aT;
                    cu[i][j][k][id_aB] = -aB;
                }

                bsrc = sC * vol + aP0 * u[i][j][k];
                cu[i][j][k][id_bsrc] = bsrc;

                // cv和cw矩阵
                cv[i][j][k][id_aP] = aP;
                cv[i][j][k][id_aE] = -aE;
                cv[i][j][k][id_aW] = -aW;
                cv[i][j][k][id_aN] = -aN;
                cv[i][j][k][id_aS] = -aS;
                if (dim == 3) {
                    cv[i][j][k][id_aT] = -aT;
                    cv[i][j][k][id_aB] = -aB;
                }

                bsrc = sC * vol + aP0 * v[i][j][k];
                cv[i][j][k][id_bsrc] = bsrc;

                if (dim == 3) {
                    cw[i][j][k][id_aP] = aP;
                    cw[i][j][k][id_aE] = -aE;
                    cw[i][j][k][id_aW] = -aW;
                    cw[i][j][k][id_aN] = -aN;
                    cw[i][j][k][id_aS] = -aS;
                    cw[i][j][k][id_aT] = -aT;
                    cw[i][j][k][id_aB] = -aB;

                    bsrc = sC * vol + aP0 * w[i][j][k];
                    cw[i][j][k][id_bsrc] = bsrc;
                }
            }
        }
    }
}


void momentum_coef_bcs(StructuredMesh& caseData, FluidBoundary& fluidBoundary, int itest, int dim, int ncx, int ncy, int ncz, int dir,
                       std::vector<std::vector<std::vector<std::vector<double>>>>& cu) {
    
    const auto& bcid = fluidBoundary.bcid;

    int fid_e = fluidBoundary.fid_e;
    int fid_w = fluidBoundary.fid_w;
    int fid_n = fluidBoundary.fid_n;
    int fid_s = fluidBoundary.fid_s;
    int fid_t = dim == 3 ? fluidBoundary.fid_t : 0;
    int fid_b = dim == 3 ? fluidBoundary.fid_b : 0;

    const auto& bcs = fluidBoundary.bcs;
    const auto& bcs_vel = fluidBoundary.bcs_vel;

    int bc_none = fluidBoundary.bc_none;
    int bc_wall = fluidBoundary.bc_wall;
    int bc_inlet = fluidBoundary.bc_inlet;
    int bc_outlet = fluidBoundary.bc_outlet;

    
    int id_aP = caseData.id_aP;
    int id_aE = caseData.id_aE;
    int id_aW = caseData.id_aW;
    int id_aN = caseData.id_aN;
    int id_aS = caseData.id_aS;
    int id_aT = caseData.id_aT;
    int id_aB = caseData.id_aB;
    int id_bsrc = caseData.id_bsrc;

    for (int k = 0; k < ncz; ++k) {
        for (int j = 0; j < ncy; ++j) {
            for (int i = 0; i < ncx; ++i) {

                int bcid_e = bcid[i][j][k][fid_e];
                int bcid_w = bcid[i][j][k][fid_w];
                int bcid_n = bcid[i][j][k][fid_n];
                int bcid_s = bcid[i][j][k][fid_s];

                int bc_e = bcs[bcid_e].type;
                int bc_w = bcs[bcid_w].type;
                int bc_n = bcs[bcid_n].type;
                int bc_s = bcs[bcid_s].type;

                double bc_ue, bc_uw, bc_un, bc_us;

                // 0,1,2 X,Y,Z 方向
                if (dir == 0) { 
                    bc_ue = bcs_vel[bcid_e].u;
                    bc_uw = bcs_vel[bcid_w].u;
                    bc_un = bcs_vel[bcid_n].u;
                    bc_us = bcs_vel[bcid_s].u;
                } else if (dir == 1) { 
                    bc_ue = bcs_vel[bcid_e].v;
                    bc_uw = bcs_vel[bcid_w].v;
                    bc_un = bcs_vel[bcid_n].v;
                    bc_us = bcs_vel[bcid_s].v;
                } else if (dir == 2) { 
                    bc_ue = bcs_vel[bcid_e].w;
                    bc_uw = bcs_vel[bcid_w].w;
                    bc_un = bcs_vel[bcid_n].w;
                    bc_us = bcs_vel[bcid_s].w;
                }

                // E
                if (bc_e == bc_inlet) {
                    cu[i][j][k][id_aP] -= 2.0 * cu[i][j][k][id_aE];
                    cu[i][j][k][id_aW] += cu[i][j][k][id_aE] / 3.0;
                    cu[i][j][k][id_bsrc] -= 8.0 * cu[i][j][k][id_aE] * bc_ue / 3.0;
                    cu[i][j][k][id_aE] = 0.0;
                } else if (bc_e == bc_outlet) {
                    cu[i][j][k][id_aP] += cu[i][j][k][id_aE];
                    cu[i][j][k][id_aE] = 0.0;
                } else if (bc_e == bc_wall) {
                    cu[i][j][k][id_aP] -= 2.0 * cu[i][j][k][id_aE];
                    cu[i][j][k][id_aW] += cu[i][j][k][id_aE] / 3.0;
                    cu[i][j][k][id_bsrc] -= 8.0 * cu[i][j][k][id_aE] * bc_ue / 3.0;
                    cu[i][j][k][id_aE] = 0.0;
                }

                // W
                if (bc_w == bc_inlet) {
                    cu[i][j][k][id_aP] -= 2.0 * cu[i][j][k][id_aW];
                    cu[i][j][k][id_aE] += cu[i][j][k][id_aW] / 3.0;
                    cu[i][j][k][id_bsrc] -= 8.0 * cu[i][j][k][id_aW] * bc_uw / 3.0;
                    cu[i][j][k][id_aW] = 0.0;
                } else if (bc_w == bc_outlet) {
                    cu[i][j][k][id_aP] += cu[i][j][k][id_aW];
                    cu[i][j][k][id_aW] = 0.0;
                } else if (bc_w == bc_wall) {
                    cu[i][j][k][id_aP] -= 2.0 * cu[i][j][k][id_aW];
                    cu[i][j][k][id_aE] += cu[i][j][k][id_aW] / 3.0;
                    cu[i][j][k][id_bsrc] -= 8.0 * cu[i][j][k][id_aW] * bc_uw / 3.0;
                    cu[i][j][k][id_aW] = 0.0;
                }

                // M
                if (bc_n == bc_inlet) {
                    cu[i][j][k][id_aP] -= 2.0 * cu[i][j][k][id_aN];
                    cu[i][j][k][id_aS] += cu[i][j][k][id_aN] / 3.0;
                    cu[i][j][k][id_bsrc] -= 8.0 * cu[i][j][k][id_aN] * bc_un / 3.0;
                    cu[i][j][k][id_aN] = 0.0;
                } else if (bc_n == bc_outlet) {
                    cu[i][j][k][id_aP] += cu[i][j][k][id_aN];
                    cu[i][j][k][id_aN] = 0.0;
                } else if (bc_n == bc_wall) {
                    cu[i][j][k][id_aP] -= 2.0 * cu[i][j][k][id_aN];
                    cu[i][j][k][id_aS] += cu[i][j][k][id_aN] / 3.0;
                    cu[i][j][k][id_bsrc] -= 8.0 * cu[i][j][k][id_aN] * bc_un / 3.0;
                    cu[i][j][k][id_aN] = 0.0;
                }

                // S
                if (bc_s == bc_inlet) {
                    cu[i][j][k][id_aP] -= 2.0 * cu[i][j][k][id_aS];
                    cu[i][j][k][id_aN] += cu[i][j][k][id_aS] / 3.0;
                    cu[i][j][k][id_bsrc] -= 8.0 * cu[i][j][k][id_aS] * bc_us / 3.0;
                    cu[i][j][k][id_aS] = 0.0;
                } else if (bc_s == bc_outlet) {
                    cu[i][j][k][id_aP] += cu[i][j][k][id_aS];
                    cu[i][j][k][id_aS] = 0.0;
                } else if (bc_s == bc_wall) {
                    cu[i][j][k][id_aP] -= 2.0 * cu[i][j][k][id_aS];
                    cu[i][j][k][id_aN] += cu[i][j][k][id_aS] / 3.0;
                    cu[i][j][k][id_bsrc] -= 8.0 * cu[i][j][k][id_aS] * bc_us / 3.0;
                    cu[i][j][k][id_aS] = 0.0;
                }

                if (dim == 3) {
                    // 3D Z方向
                }
            }
        }
    }
}


void momentum_gradp(StructuredMesh& caseData, FluidBoundary& fluidBoundary, int itest, int dim, int ncx, int ncy, int ncz,
                    const std::vector<std::vector<std::vector<double>>>& p,
                    std::vector<std::vector<std::vector<std::vector<double>>>>& cu,
                    std::vector<std::vector<std::vector<std::vector<double>>>>& cv,
                    std::vector<std::vector<std::vector<std::vector<double>>>>& cw) {

    const auto& bcid = fluidBoundary.bcid;

    int fid_e = fluidBoundary.fid_e;
    int fid_w = fluidBoundary.fid_w;
    int fid_n = fluidBoundary.fid_n;
    int fid_s = fluidBoundary.fid_s;
    int fid_t = dim == 3 ? fluidBoundary.fid_t : 0;
    int fid_b = dim == 3 ? fluidBoundary.fid_b : 0;

    const auto& bcs = fluidBoundary.bcs;

    int bc_none = fluidBoundary.bc_none;
    int bc_wall = fluidBoundary.bc_wall;
    int bc_inlet = fluidBoundary.bc_inlet;
    int bc_outlet = fluidBoundary.bc_outlet;

    double p_outlet = fluidBoundary.p_outlet;


    double dx = caseData.dx;
    double dy = caseData.dy;
    double dz = caseData.dz;

    double area_x = dy * dz;
    double area_y = dx * dz;
    double area_z = dx * dy;

    for (int k = 0; k < ncz; ++k) {
        for (int j = 0; j < ncy; ++j) {
            for (int i = 0; i < ncx; ++i) {

                int bcid_e = bcid[i][j][k][fid_e];
                int bcid_w = bcid[i][j][k][fid_w];
                int bcid_n = bcid[i][j][k][fid_n];
                int bcid_s = bcid[i][j][k][fid_s];

                int bc_e = bcs[bcid_e].type;
                int bc_w = bcs[bcid_w].type;
                int bc_n = bcs[bcid_n].type;
                int bc_s = bcs[bcid_s].type;

                double pl, pr;

                // x方向压力梯度
                if (bc_e == bc_none && bc_w == bc_none) { 
                    if (i > 0 && i < ncx - 1) {
                        pl = p[i - 1][j][k];
                        pr = p[i + 1][j][k];
                    } else {
                        pl = p[i][j][k];
                        pr = p[i][j][k];
                    }
                } else {
                    if (bc_e == bc_inlet || bc_e == bc_wall) {
                        pl = p[i - 1][j][k];
                        pr = p[i][j][k];
                    } else if (bc_e == bc_outlet) {
                        pl = p[i - 1][j][k];
                        pr = p_outlet;
                    }

                    if (bc_w == bc_inlet || bc_w == bc_wall) {
                        pl = p[i][j][k];
                        pr = p[i + 1][j][k];
                    } else if (bc_w == bc_outlet) {
                        pl = p_outlet;
                        pr = p[i + 1][j][k];
                    }
                }

                cu[i][j][k][caseData.id_bsrc] += 0.5 * (pl - pr) * area_x;

                // y
                if (bc_n == bc_none && bc_s == bc_none) {  
                    if (j > 0 && j < ncy - 1) {
                        pl = p[i][j - 1][k];
                        pr = p[i][j + 1][k];
                    } else {
                        pl = p[i][j][k];
                        pr = p[i][j][k];
                    }
                } else {
                    if (bc_n == bc_inlet || bc_n == bc_wall) {
                        pl = p[i][j - 1][k];
                        pr = p[i][j][k];
                    } else if (bc_n == bc_outlet) {
                        pl = p[i][j - 1][k];
                        pr = p_outlet;
                    }

                    if (bc_s == bc_inlet || bc_s == bc_wall) {
                        pl = p[i][j][k];
                        pr = p[i][j + 1][k];
                    } else if (bc_s == bc_outlet) {
                        pl = p_outlet;
                        pr = p[i][j + 1][k];
                    }
                }

                cv[i][j][k][caseData.id_bsrc] += 0.5 * (pl - pr) * area_y;

                // z方向的压力梯度
                if (dim == 3) {

                }
            }
        }
    }
}


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
                             const std::vector<std::vector<std::vector<std::vector<double>>>>& cw) {

    const auto& bcid = fluidBoundary.bcid;

    int fid_e = fluidBoundary.fid_e;
    int fid_w = fluidBoundary.fid_w;
    int fid_n = fluidBoundary.fid_n;
    int fid_s = fluidBoundary.fid_s;
    int fid_t = dim == 3 ? fluidBoundary.fid_t : 0;
    int fid_b = dim == 3 ? fluidBoundary.fid_b : 0;

    const auto& bcs = fluidBoundary.bcs;

    int bc_none = fluidBoundary.bc_none;
    int bc_wall = fluidBoundary.bc_wall;
    int bc_inlet = fluidBoundary.bc_inlet;
    int bc_outlet = fluidBoundary.bc_outlet;

    double p_outlet = fluidBoundary.p_outlet;


    double dx = caseData.dx;
    double dy = caseData.dy;
    double dz = caseData.dz;

    double area_x = dy * dz;
    double area_y = dx * dz;
    double area_z = dx * dy;
    double vol = dx * dy * dz;

    double idt = 1.0 / dt;

    // x方向的面速度修正
    for (int k = 0; k < ncz; ++k) {
        for (int j = 0; j < ncy; ++j) {
            for (int i = 1; i < ncx; ++i) {

                int bcid_e = bcid[i][j][k][fid_e];
                int bcid_w = bcid[i - 1][j][k][fid_w];
                int bc_e = bcs[bcid_e].type;
                int bc_w = bcs[bcid_w].type;

                double uf0 = uf[i][j][k];

                double pe = p[i][j][k];
                double pw = p[i - 1][j][k];

                double pww, pee;

                if (i == 1) {
                    pww = 0.0;
                } else {
                    pww = p[i - 2][j][k];
                }

                if (i == ncx - 1) {
                    pee = 0.0;
                } else {
                    pee = p[i + 1][j][k];
                }

                if (bc_w == bc_inlet || bc_w == bc_wall) {
                    pww = p[i - 1][j][k];
                } else if (bc_w == bc_outlet) {
                    pww = p_outlet;
                }

                if (bc_e == bc_inlet || bc_e == bc_wall) {
                    pee = p[i][j][k];
                } else if (bc_e == bc_outlet) {
                    pee = p_outlet;
                }

                double dpm = 0.5 * (pe - pww);
                double dpp = 0.5 * (pee - pw);
                double dpf = pe - pw;

                double uf_new = (
                    (cu[i - 1][j][k][caseData.id_aP] * u[i - 1][j][k] + cu[i][j][k][caseData.id_aP] * u[i][j][k])
                    / (cu[i - 1][j][k][caseData.id_aP] + cu[i][j][k][caseData.id_aP])
                    + (
                        0.5 * dpm / cu[i - 1][j][k][caseData.id_aP]
                        + 0.5 * dpp / cu[i][j][k][caseData.id_aP]
                        - 0.5 * dpf * (1.0 / cu[i - 1][j][k][caseData.id_aP] + 1.0 / cu[i][j][k][caseData.id_aP])
                    ) * area_x
                );

                double cf0 = dens[i][j][k] * idt;
                double cf = dens[i][j][k] * idt;
                double df = 0.5 * (vol / cu[i - 1][j][k][caseData.id_aP] + vol / cu[i][j][k][caseData.id_aP]);
                double df_hat = df / (1.0 + cf * df);
                uf[i][j][k] = uf_new + cf0 * df_hat * (uf0 - 0.5 * (u0[i - 1][j][k] + u0[i][j][k]));
            }
        }
    }

    // y方向
    for (int k = 0; k < ncz; ++k) {
        for (int j = 1; j < ncy; ++j) {
            for (int i = 0; i < ncx; ++i) {

                int bcid_n = bcid[i][j][k][fid_n];
                int bcid_s = bcid[i][j - 1][k][fid_s];
                int bc_n = bcs[bcid_n].type;
                int bc_s = bcs[bcid_s].type;

                double vf0 = vf[i][j][k];

                double pe = p[i][j][k];
                double pw = p[i][j - 1][k];

                double pww, pee;

                if (j == 1) {
                    pww = 0.0;
                } else {
                    pww = p[i][j - 2][k];
                }

                if (j == ncy - 1) {
                    pee = 0.0;
                } else {
                    pee = p[i][j + 1][k];
                }

                if (bc_s == bc_inlet || bc_s == bc_wall) {
                    pww = p[i][j - 1][k];
                } else if (bc_s == bc_outlet) {
                    pww = p_outlet;
                }

                if (bc_n == bc_inlet || bc_n == bc_wall) {
                    pee = p[i][j][k];
                } else if (bc_n == bc_outlet) {
                    pee = p_outlet;
                }

                double dpm = 0.5 * (pe - pww);
                double dpp = 0.5 * (pee - pw);
                double dpf = pe - pw;

                double vf_new = (
                    (cv[i][j - 1][k][caseData.id_aP] * v[i][j - 1][k] + cv[i][j][k][caseData.id_aP] * v[i][j][k])
                    / (cv[i][j - 1][k][caseData.id_aP] + cv[i][j][k][caseData.id_aP])
                    + (
                        0.5 * dpm / cv[i][j - 1][k][caseData.id_aP]
                        + 0.5 * dpp / cv[i][j][k][caseData.id_aP]
                        - 0.5 * dpf * (1.0 / cv[i][j - 1][k][caseData.id_aP] + 1.0 / cv[i][j][k][caseData.id_aP])
                    ) * area_y
                );

                double cf0 = dens[i][j][k] * idt;
                double cf = dens[i][j][k] * idt;
                double df = 0.5 * (vol / cv[i][j - 1][k][caseData.id_aP] + vol / cv[i][j][k][caseData.id_aP]);
                double df_hat = df / (1.0 + cf * df);
                vf[i][j][k] = vf_new + cf0 * df_hat * (vf0 - 0.5 * (v0[i][j - 1][k] + v0[i][j][k]));
            }
        }
    }

    // z 方向
    if (dim == 3) {

    }
}


void pressure_coefs(StructuredMesh& caseData, FluidBoundary& fluidBoundary, int itest, int dim, int ncx, int ncy, int ncz,
                    const std::vector<std::vector<std::vector<double>>>& dens,
                    const std::vector<std::vector<std::vector<double>>>& uf,
                    const std::vector<std::vector<std::vector<double>>>& vf,
                    const std::vector<std::vector<std::vector<double>>>& wf,
                    const std::vector<std::vector<std::vector<std::vector<double>>>>& cu,
                    const std::vector<std::vector<std::vector<std::vector<double>>>>& cv,
                    const std::vector<std::vector<std::vector<std::vector<double>>>>& cw,
                    std::vector<std::vector<std::vector<std::vector<double>>>>& cp) {

    const auto& bcid = fluidBoundary.bcid;

    int fid_e = fluidBoundary.fid_e;
    int fid_w = fluidBoundary.fid_w;
    int fid_n = fluidBoundary.fid_n;
    int fid_s = fluidBoundary.fid_s;
    int fid_t = dim == 3 ? fluidBoundary.fid_t : 0;
    int fid_b = dim == 3 ? fluidBoundary.fid_b : 0;

    const auto& bcs = fluidBoundary.bcs;

    int bc_none = fluidBoundary.bc_none;
    int bc_wall = fluidBoundary.bc_wall;
    int bc_inlet = fluidBoundary.bc_inlet;
    int bc_outlet = fluidBoundary.bc_outlet;

    double p_outlet = fluidBoundary.p_outlet;


    double dx = caseData.dx;
    double dy = caseData.dy;
    double dz = caseData.dz;

    double area_x = dy * dz;
    double area_y = dx * dz;
    double area_z = dx * dy;
    double vol = dx * dy * dz;

    for (int k = 0; k < ncz; ++k) {
        for (int j = 0; j < ncy; ++j) {
            for (int i = 0; i < ncx; ++i) {
                
                int bcid_e = bcid[i][j][k][fid_e];
                int bcid_w = bcid[i][j][k][fid_w];
                int bcid_n = bcid[i][j][k][fid_n];
                int bcid_s = bcid[i][j][k][fid_s];

                int bc_e = bcs[bcid_e].type;
                int bc_w = bcs[bcid_w].type;
                int bc_n = bcs[bcid_n].type;
                int bc_s = bcs[bcid_s].type;

                // 系数初始化
                double aE = 0.0, aW = 0.0;
                double aN = 0.0, aS = 0.0;
                double aT = 0.0, aB = 0.0;
                double aP = 0.0, bsrc = 0.0;

                double rho_w = 0.0, rho_e = 0.0;
                double rho_s = 0.0, rho_n = 0.0;

                // 计算东侧面通量
                if (bc_e != bc_none || i == ncx - 1) {
                    aE = 0.0;
                    rho_e = dens[i][j][k];
                } else {
                    rho_e = 0.5 * (dens[i][j][k] + dens[i + 1][j][k]);
                    aE = rho_e * area_x * area_x * 0.5 * (1.0 / cu[i + 1][j][k][caseData.id_aP] + 1.0 / cu[i][j][k][caseData.id_aP]);
                }

                // W
                if (bc_w != bc_none || i == 0) {
                    aW = 0.0;
                    rho_w = dens[i][j][k];
                } else {
                    rho_w = 0.5 * (dens[i][j][k] + dens[i - 1][j][k]);
                    aW = rho_w * area_x * area_x * 0.5 * (1.0 / cu[i - 1][j][k][caseData.id_aP] + 1.0 / cu[i][j][k][caseData.id_aP]);
                }


                // S
                if (bc_s != bc_none || j == 0) {
                    aS = 0.0;
                    rho_s = dens[i][j][k];
                } else {
                    rho_s = 0.5 * (dens[i][j][k] + dens[i][j - 1][k]);
                    aS = rho_s * area_y * area_y * 0.5 * (1.0 / cv[i][j - 1][k][caseData.id_aP] + 1.0 / cv[i][j][k][caseData.id_aP]);
                }

                // N
                if (bc_n != bc_none || j == ncy - 1) {
                    aN = 0.0;
                    rho_n = dens[i][j][k];
                } else {
                    rho_n = 0.5 * (dens[i][j][k] + dens[i][j + 1][k]);
                    aN = rho_n * area_y * area_y * 0.5 * (1.0 / cv[i][j + 1][k][caseData.id_aP] + 1.0 / cv[i][j][k][caseData.id_aP]);
                }

                // 顶部和底部
                if (dim == 3) {

                }

                aP = aE + aW + aN + aS + aT + aB;

                bsrc = (rho_w * uf[i][j][k] - rho_e * uf[i + 1][j][k]) * area_x +
                       (rho_s * vf[i][j][k] - rho_n * vf[i][j + 1][k]) * area_y;

                cp[i][j][k][caseData.id_aP] = aP;
                cp[i][j][k][caseData.id_aE] = -aE;
                cp[i][j][k][caseData.id_aW] = -aW;
                cp[i][j][k][caseData.id_aN] = -aN;
                cp[i][j][k][caseData.id_aS] = -aS;

                if (dim == 3) {
                    cp[i][j][k][caseData.id_aT] = -aT;
                    cp[i][j][k][caseData.id_aB] = -aB;
                }

                cp[i][j][k][caseData.id_bsrc] = bsrc;
            }
        }
    }
}


void correct_pressure(int ncx, int ncy, int ncz, double relax_p,
                      std::vector<std::vector<std::vector<double>>>& pp,
                      std::vector<std::vector<std::vector<double>>>& p) {
    for (int k = 0; k < ncz; ++k) {
        for (int j = 0; j < ncy; ++j) {
            for (int i = 0; i < ncx; ++i) {
                p[i][j][k] += relax_p * pp[i][j][k];
            }
        }
    }
}


void correct_velocity(int dim, int ncx, int ncy, int ncz,
                      const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z,
                      const std::vector<std::vector<std::vector<std::vector<double>>>>& cu,
                      const std::vector<std::vector<std::vector<std::vector<double>>>>& cv,
                      const std::vector<std::vector<std::vector<double>>>& pp,
                      std::vector<std::vector<std::vector<double>>>& u,
                      std::vector<std::vector<std::vector<double>>>& v,
                      std::vector<std::vector<std::vector<double>>>& uf,
                      std::vector<std::vector<std::vector<double>>>& vf) {
    double dx = x[1] - x[0];
    double dy = y[1] - y[0];
    double dz = (dim == 3) ? (z[1] - z[0]) : 1.0;

    double area_x = dy * dz;
    double area_y = dx * dz;

    // 修正u速度
    for (int k = 0; k < ncz; ++k) {
        for (int j = 0; j < ncy; ++j) {
            for (int i = 0; i < ncx; ++i) {
                double pl = (i == 0) ? 0.0 : pp[i - 1][j][k];
                double pr = (i == ncx - 1) ? 0.0 : pp[i + 1][j][k];

                double d = area_x * 0.5 / cu[i][j][k][0];
                u[i][j][k] += d * (pl - pr);
            }
        }
    }

    // 修正v
    for (int k = 0; k < ncz; ++k) {
        for (int j = 0; j < ncy; ++j) {
            for (int i = 0; i < ncx; ++i) {
                double pl = (j == 0) ? 0.0 : pp[i][j - 1][k];
                double pr = (j == ncy - 1) ? 0.0 : pp[i][j + 1][k];

                double d = area_y * 0.5 / cv[i][j][k][0];
                v[i][j][k] += d * (pl - pr);
            }
        }
    }

    // 修正uf
    for (int k = 0; k < ncz; ++k) {
        for (int j = 0; j < ncy; ++j) {
            for (int i = 1; i < ncx; ++i) {
                double d = area_x * 0.5 * (1.0 / cu[i - 1][j][k][0] + 1.0 / cu[i][j][k][0]);
                uf[i][j][k] += d * (pp[i - 1][j][k] - pp[i][j][k]);
            }
        }
    }

    // 修正vf
    for (int k = 0; k < ncz; ++k) {
        for (int j = 1; j < ncy; ++j) {
            for (int i = 0; i < ncx; ++i) {
                double d = area_y * 0.5 * (1.0 / cv[i][j - 1][k][0] + 1.0 / cv[i][j][k][0]);
                vf[i][j][k] += d * (pp[i][j - 1][k] - pp[i][j][k]);
            }
        }
    }
}

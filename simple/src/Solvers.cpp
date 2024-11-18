#include "Solvers.h"


std::pair<double, double> eqn_scalar_norm2(StructuredMesh& caseData, int dim, int it_nl, int ncx, int ncy, int ncz,
                                           const std::vector<std::vector<std::vector<double>>>& u0,
                                           const std::vector<std::vector<std::vector<double>>>& u,
                                           const std::string& var) {

    double l2_u = 0.0;
    double l2_max_u = 0.0;

    // 根据不同的变量名存入目标变量
    if (var == "temp") {
        l2_u = caseData.l2_t;
        l2_max_u = caseData.l2_max_t;
    } else if (var == "p") {
        l2_u = caseData.l2_p;
        l2_max_u = caseData.l2_max_p;
    } else if (var == "u") {
        l2_u = caseData.l2_u;
        l2_max_u = caseData.l2_max_u;
    } else if (var == "v") {
        l2_u = caseData.l2_v;
        l2_max_u = caseData.l2_max_v;
    } else if (var == "w") {
        l2_u = caseData.l2_w;
        l2_max_u = caseData.l2_max_w;
    } else {
        std::cerr << "Unknown variable: " << var << std::endl;
        return {0.0, 0.0};
    }

    // 计算 u 和 u0 之间的平方差和
    double sum_sq_diff = 0.0;
    int ncell = ncx * ncy * ncz;                 // 总网格单元数
    for (int k = 0; k < ncz; ++k) {
        for (int j = 0; j < ncy; ++j) {
            for (int i = 0; i < ncx; ++i) {
                double diff = u[i][j][k] - u0[i][j][k];
                sum_sq_diff += diff * diff;
            }
        }
    }

    // L2范数计算
    l2_u = std::sqrt(sum_sq_diff / static_cast<double>(ncell));
    l2_max_u = std::max(l2_u, l2_max_u);

    // 数据更新
    if (var == "temp") {
        caseData.l2_t = l2_u;
        caseData.l2_max_t = l2_max_u;
    } else if (var == "p") {
        caseData.l2_p = l2_u;
        caseData.l2_max_p = l2_max_u;
    } else if (var == "u") {
        caseData.l2_u = l2_u;
        caseData.l2_max_u = l2_max_u;
    } else if (var == "v") {
        caseData.l2_v = l2_u;
        caseData.l2_max_v = l2_max_u;
    } else if (var == "w") {
        caseData.l2_w = l2_u;
        caseData.l2_max_w = l2_max_u;
    }

    return {l2_u, l2_max_u};
}


void eqn_scalar_norm2_cp(StructuredMesh& caseData, int dim, int it_nl, int ncx, int ncy, int ncz,
                         const std::vector<std::vector<std::vector<std::vector<double>>>>& c,
                         const std::string& var) {

    int id_bsrc = caseData.id_bsrc;
    double l2_u = 0.0;
    double l2_max_u = 0.0;

    if (var == "cp") {
        l2_u = caseData.l2_pp;
        l2_max_u = caseData.l2_max_pp;
    } else {
        std::cerr << "Unknown variable: " << var << std::endl;
        return;
    }

    
    double sum_c4 = 0.0;
    int ncell = ncx * ncy * ncz;
    for (int k = 0; k < ncz; ++k) {
        for (int j = 0; j < ncy; ++j) {
            for (int i = 0; i < ncx; ++i) {
                double c_val = c[i][j][k][id_bsrc];
                double c_sq = c_val * c_val;
                sum_c4 += c_sq * c_sq;
            }
        }
    }

    l2_u = std::sqrt(sum_c4 / static_cast<double>(ncell));
    l2_max_u = std::max(l2_u, l2_max_u);

    if (var == "cp") {
        caseData.l2_pp = l2_u;
        caseData.l2_max_pp = l2_max_u;
    }
}


void scalar_pj(StructuredMesh& caseData, int dim, int it_nl, int niter, double relax,
               int ncx, int ncy, int ncz, int ncoef,
               const std::vector<std::vector<std::vector<std::vector<double>>>>& ct,
               std::vector<std::vector<std::vector<double>>>& t,
               bool initzero, double res) {

    // 提取系数矩阵的索引
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

    int nsteps = caseData.nsteps;

    int res_freq = 1;

    // 设置初始解为零
    if (initzero) {
        for (int k = 0; k < ncz; ++k)
            for (int j = 0; j < ncy; ++j)
                for (int i = 0; i < ncx; ++i)
                    t[i][j][k] = 0.0;
    }

    std::vector<std::vector<std::vector<double>>> t0 = t;

    double max_norm = -1e20;        // 最大
    double initial_norm = 0.0;      // 初始

    for (int it = 1; it <= niter; ++it) {

        t0 = t;

        double norm2 = 0.0;

        for (int k = 0; k < ncz; ++k) {
            for (int j = 0; j < ncy; ++j) {
                for (int i = 0; i < ncx; ++i) {
                     // 获取相邻网格的值
                    double tw = (i == 0) ? 0.0 : t0[i - 1][j][k];
                    double te = (i == ncx - 1) ? 0.0 : t0[i + 1][j][k];
                    double ts = (j == 0) ? 0.0 : t0[i][j - 1][k];
                    double tn = (j == ncy - 1) ? 0.0 : t0[i][j + 1][k];

                    double tt = 0.0, tb = 0.0;
                    if (dim == 3) {
                        tb = (k == 0) ? 0.0 : t0[i][j][k - 1];
                        tt = (k == ncz - 1) ? 0.0 : t0[i][j][k + 1];
                    }

                    // 计算新的值
                    double t_new = (
                        - ct[i][j][k][id_aE] * te
                        - ct[i][j][k][id_aW] * tw
                        - ct[i][j][k][id_aN] * tn
                        - ct[i][j][k][id_aS] * ts
                        + ct[i][j][k][id_bsrc]
                    );

                    if (dim == 3) {
                        t_new -= ct[i][j][k][id_aT] * tt + ct[i][j][k][id_aB] * tb;
                    }

                    t_new = t_new / ct[i][j][k][id_aP];

                    // 更新松弛因子
                    double du = relax * (t_new - t0[i][j][k]);
                    t[i][j][k] = t0[i][j][k] + du;
                    norm2 += du * du;
                }
            }
        }

        // 计算当前范数
        int ncell = ncx * ncy * ncz;
        norm2 = std::sqrt(norm2 / static_cast<double>(ncell));
        if (it == 1) {
            initial_norm = norm2 + 1e-20;
        }
        max_norm = std::max(norm2, max_norm) + 1e-20;

        // 判断收敛条件
        double rel_norm = norm2 / max_norm;
        if (rel_norm < res || it == niter) {
            caseData.total_linsol_iters += it;
            if (it_nl % res_freq == 0 || it_nl == 1 || it_nl == nsteps) {
                std::cout << "it_nl, it, tot_it, norm2, init, max, rel "
                          << it_nl << " " << it << " " << caseData.total_linsol_iters << " "
                          << norm2 << " " << initial_norm << " " << max_norm << " " << rel_norm << std::endl;
            }
            break;
        }
    }
}


void scalar_gs(StructuredMesh& caseData, int dim, int it_nl, int niter, double relax,
               int ncx, int ncy, int ncz, int ncoef,
               const std::vector<std::vector<std::vector<std::vector<double>>>>& ct,
               std::vector<std::vector<std::vector<double>>>& t,
               bool initzero, double res) {

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

    int nsteps = caseData.nsteps;

    int res_freq = 1;

    if (initzero) {
        for (int k = 0; k < ncz; ++k)
            for (int j = 0; j < ncy; ++j)
                for (int i = 0; i < ncx; ++i)
                    t[i][j][k] = 0.0;
    }

    double max_norm = -1e20;
    double initial_norm = 0.0;

    for (int it = 1; it <= niter; ++it) {

        double norm2 = 0.0;

        for (int k = 0; k < ncz; ++k) {
            for (int j = 0; j < ncy; ++j) {
                for (int i = 0; i < ncx; ++i) {
                    double tw = (i == 0) ? 0.0 : t[i - 1][j][k];
                    double te = (i == ncx - 1) ? 0.0 : t[i + 1][j][k];
                    double ts = (j == 0) ? 0.0 : t[i][j - 1][k];
                    double tn = (j == ncy - 1) ? 0.0 : t[i][j + 1][k];

                    double tt = 0.0, tb = 0.0;
                    if (dim == 3) {
                        tb = (k == 0) ? 0.0 : t[i][j][k - 1];
                        tt = (k == ncz - 1) ? 0.0 : t[i][j][k + 1];
                    }

                    double t_new = (
                        - ct[i][j][k][id_aE] * te
                        - ct[i][j][k][id_aW] * tw
                        - ct[i][j][k][id_aN] * tn
                        - ct[i][j][k][id_aS] * ts
                        + ct[i][j][k][id_bsrc]
                    );

                    if (dim == 3) {
                        t_new -= ct[i][j][k][id_aT] * tt + ct[i][j][k][id_aB] * tb;
                    }

                    t_new = t_new / ct[i][j][k][id_aP];

                    double du = relax * (t_new - t[i][j][k]);
                    t[i][j][k] += du;
                    norm2 += du * du;
                }
            }
        }

        int ncell = ncx * ncy * ncz;
        norm2 = std::sqrt(norm2 / static_cast<double>(ncell));

        if (it == 1) {
            initial_norm = norm2 + 1e-20;
        }

        max_norm = std::max(norm2, max_norm) + 1e-20;

        double rel_norm = norm2 / max_norm;

        if (rel_norm < res || it == niter) {
            caseData.total_linsol_iters += it;
            if (it_nl % res_freq == 0 || it_nl == 1 || it_nl == nsteps) {
                std::cout << "it_nl, it, tot_it, norm2, init, max, rel "
                          << it_nl << " " << it << " " << caseData.total_linsol_iters << " "
                          << norm2 << " " << initial_norm << " " << max_norm << " " << rel_norm << std::endl;
            }
            break;
        }
    }
}

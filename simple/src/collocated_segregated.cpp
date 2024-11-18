#include "collocated_segregated.h"

bool collocated_segregated(int it, StructuredMesh& caseData, Fluid& fluid, FluidBoundary& fluidBoundary) {
    
    // 获取维度
    int dim = caseData.dim;

    // 更新解
    caseData.t0 = caseData.t;
    caseData.u0 = caseData.u;
    caseData.v0 = caseData.v;
    caseData.w0 = caseData.w;
    caseData.p0 = caseData.p;

    // 网格信息
    int ncx = caseData.ncx;
    int ncy = caseData.ncy;
    int ncz = caseData.ncz;

    // 系数矩阵
    int ncoef = caseData.ncoef;
    int ncoef_p = caseData.ncoef_p;

    // 时间步长
    double dt = caseData.dt;

    // 变量
    auto& x = caseData.x;
    auto& y = caseData.y;
    auto& z = caseData.z;

    auto& t = caseData.t;
    auto& t0 = caseData.t0;
    auto& u = caseData.u;
    auto& u0 = caseData.u0;
    auto& v = caseData.v;
    auto& v0 = caseData.v0;
    auto& w = caseData.w;
    auto& w0 = caseData.w0;
    auto& p = caseData.p;
    auto& p0 = caseData.p0;
    auto& pp = caseData.pp;

    auto& uf = caseData.uf;
    auto& vf = caseData.vf;
    auto& wf = caseData.wf;

    auto& ct = caseData.ct;
    auto& cu = caseData.cu;
    auto& cv = caseData.cv;
    auto& cw = caseData.cw;
    auto& cp = caseData.cp;

    // 方程标识与使用
    int ieqn = caseData.ieqn;
    int eqn_conduction = caseData.eqn_conduction;
    int eqn_flow = caseData.eqn_flow;
    int eqn_conduction_flow = caseData.eqn_conduction_flow;

    //求解方法选择
    int isolver = caseData.isolver;
    int solver_pj = caseData.solver_pj;
    int solver_gs = caseData.solver_gs;

    // 迭代与松弛
    int niter_t = caseData.niter_t;
    double relax_t = caseData.relax_t;
    int niter_u = caseData.niter_u;
    double relax_u = caseData.relax_u;
    int niter_v = caseData.niter_v;
    double relax_v = caseData.relax_v;
    int niter_w = caseData.niter_w;
    double relax_w = caseData.relax_w;
    int niter_p = caseData.niter_p;
    double relax_p = caseData.relax_p;
    double res_t = caseData.res_t;
    double res_u = caseData.res_u;
    double res_v = caseData.res_v;
    double res_w = caseData.res_w;
    double res_p = caseData.res_p;

    if (dim == 3) {
        int niter_w = caseData.niter_w;
        double relax_w = caseData.relax_w;
        double res_w = caseData.res_w;
    }

    // tolerence
    double temp_tol = caseData.temp_tol;
    double mass_tol = caseData.mass_tol;
    bool stop_sim = caseData.stop_sim;

    // 流体热属性
    double spht = fluid.spht;
    double con = fluid.con;
    auto& dens = fluid.dens;
    auto& mu = fluid.mu;
    double heat_src = fluid.heat_src;

    //  求解类型为 流场 或者 流场、温度场同时求解（实际上这段代码只解流场）
    if (ieqn == eqn_flow || ieqn == eqn_conduction_flow) {
        set_face_vel_bc(caseData, fluidBoundary, caseData.itest, dim, ncx, ncy, ncz, u, v, w, uf, vf, wf);

        momentum_coefs(caseData, dim, ncx, ncy, ncz, ncoef, dt, mu, dens, u, v, w, uf, vf, wf, cu, cv, cw);

        for (int dir = 0; dir < dim; dir++) {
            momentum_coef_bcs(caseData, fluidBoundary, caseData.itest, dim, ncx, ncy, ncz, dir, dir == 0 ? cu : dir == 1 ? cv : cw);
        }

        momentum_gradp(caseData, fluidBoundary, caseData.itest, dim, ncx, ncy, ncz, p, cu, cv, cw);

        if (isolver == solver_pj) {
            scalar_pj(caseData, dim, it, niter_u, relax_u, ncx, ncy, ncz, ncoef, cu, u, false, res_u);
            scalar_pj(caseData, dim, it, niter_v, relax_v, ncx, ncy, ncz, ncoef, cv, v, false, res_v);
            if (dim == 3) {
                scalar_pj(caseData, dim, it, niter_w, relax_w, ncx, ncy, ncz, ncoef, cw, w, false, res_w);
            }
        } else if (isolver == solver_gs) {
            scalar_gs(caseData, dim, it, niter_u, relax_u, ncx, ncy, ncz, ncoef, cu, u, false, res_u);
            scalar_gs(caseData, dim, it, niter_v, relax_v, ncx, ncy, ncz, ncoef, cv, v, false, res_v);
            if (dim == 3) {
                scalar_gs(caseData, dim, it, niter_w, relax_w, ncx, ncy, ncz, ncoef, cw, w, false, res_w);
            }
        }

        rhie_chow_face_velocity(caseData, fluidBoundary, caseData.itest, dim, ncx, ncy, ncz, ncoef, x, y, z, dt, dens, p, u, v, w, u0, v0, w0, uf, vf, wf, cu, cv, cw);

        pressure_coefs(caseData, fluidBoundary, caseData.itest, dim, ncx, ncy, ncz, dens, uf, vf, wf, cu, cv, cw, cp);

        if (isolver == solver_pj) {
            scalar_pj(caseData, dim, it, niter_p, relax_p, ncx, ncy, ncz, ncoef_p, cp, pp, true, res_p);
        } else if (isolver == solver_gs) {
            scalar_gs(caseData, dim, it, niter_p, relax_p, ncx, ncy, ncz, ncoef_p, cp, pp, true, res_p);
        }

        correct_pressure(ncx, ncy, ncz, relax_p, caseData.pp, caseData.p);

        correct_velocity(dim, ncx, ncy, ncz, caseData.x, caseData.y, caseData.z,
                 caseData.cu, caseData.cv, caseData.pp,
                 caseData.u, caseData.v, caseData.uf, caseData.vf);

        auto [l2_u, l2_max_u] = eqn_scalar_norm2(caseData, dim, it, ncx, ncy, ncz, u0, u, "u");

        stop_sim = (l2_u / mass_tol) < mass_tol;
    }

    return stop_sim;
}

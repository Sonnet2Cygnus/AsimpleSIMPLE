// StructuredMesh.cpp

#include "StructuredMesh.h"

// 初始化
StructuredMesh::StructuredMesh()
    : dim(0), ncx(0), ncy(0), ncz(0), nx(0), ny(0), nz(0),
      dx(0.0), dy(0.0), dz(0.0),
      xmin(0.0), xmax(0.0), ymin(0.0), ymax(0.0), zmin(0.0), zmax(0.0),
      initial_uf(0.0), initial_vf(0.0), initial_wf(0.0),
      ncoef(0), ncoef_p(0),
      stop_sim(false), total_linsol_iters(0),
      l2_curr(0.0), l2_max(-1e20), l2_max_u(-1e20), l2_max_v(-1e20), l2_max_w(-1e20),
      l2_max_p(-1e20), l2_max_pp(-1e20), l2_max_t(-1e20),
      l2_u(0.0), l2_v(0.0), l2_w(0.0), l2_p(0.0), l2_pp(0.0), l2_t(0.0){

        // 初始化重力向量
        grav.clear();
}


// 创建网格
void StructuredMesh::CreateMesh(int dim, int ncx, int ncy, int ncz) {

    this->dim = dim;
    this->ncx = ncx;
    this->ncy = ncy;
    this->ncz = ncz;

    if (this->dim == 2) {
        this->ncz = 1;
    }

    // 设置网格节点数 +1
    this->nx = this->ncx + 1;
    this->ny = this->ncy + 1;
    this->nz = this->ncz + 1;

    // 输出
    std::cout << "nx, ny, nz = " << nx << ", " << ny << ", " << nz << std::endl;
    std::cout << "ncx, ncy, ncz = " << ncx << ", " << ncy << ", " << ncz << std::endl;
}

// 创建坐标
void StructuredMesh::CreateCoordinates(double xmin, double xmax, double ymin, double ymax,double zmin, double zmax) {

    // 设置边界
    this->xmin = xmin;
    this->xmax = xmax;
    this->ymin = ymin;
    this->ymax = ymax;
    this->zmin = zmin;
    this->zmax = zmax;

    if (this->dim == 2) {
        this->zmin = 0.0;
        this->zmax = 1.0;
    }

    // 初始化网格节点和中心点坐标向量
    x.resize(nx, 0.0);
    y.resize(ny, 0.0);
    z.resize(nz, 0.0);

    xc.resize(ncx, 0.0);
    yc.resize(ncy, 0.0);
    zc.resize(ncz, 0.0);

    // 计算 各方向节点坐标及单元尺寸
    dx = (xmax - xmin) / static_cast<double>(nx - 1);
    for (int i = 0; i < nx; ++i) {
        x[i] = xmin + static_cast<double>(i) * dx;
    }

    dy = (ymax - ymin) / static_cast<double>(ny - 1);
    for (int i = 0; i < ny; ++i) {
        y[i] = ymin + static_cast<double>(i) * dy;
    }

    dz = (zmax - zmin) / static_cast<double>(nz - 1);
    for (int i = 0; i < nz; ++i) {
        z[i] = zmin + static_cast<double>(i) * dz;
    }

    // 计算单元中心坐标
    for (int i = 0; i < ncx; ++i) {
        xc[i] = 0.5 * (x[i] + x[i + 1]);
    }
    for (int i = 0; i < ncy; ++i) {
        yc[i] = 0.5 * (y[i] + y[i + 1]);
    }
    for (int i = 0; i < ncz; ++i) {
        zc[i] = 0.5 * (z[i] + z[i + 1]);
    }

    // 输出
    std::cout << "dx, dy, dz = " << dx << ", " << dy << ", " << dz << std::endl;
    std::cout << "bbox xmin, xmax = " << xmin << ", " << xmax << std::endl;
    std::cout << "bbox ymin, ymax = " << ymin << ", " << ymax << std::endl;
    std::cout << "bbox zmin, zmax = " << zmin << ", " << zmax << std::endl;
    std::cout << "bbox Lx, Ly, Lz = " << xmax - xmin << ", " << ymax - ymin << ", " << zmax - zmin << std::endl;
}

// 创建场数据
void StructuredMesh::CreateFieldMeshData() {

    // u, v, w, p, t
    u.resize(ncx, std::vector<std::vector<double>>(ncy, std::vector<double>(ncz, 0.0)));
    v.resize(ncx, std::vector<std::vector<double>>(ncy, std::vector<double>(ncz, 0.0)));
    w.resize(ncx, std::vector<std::vector<double>>(ncy, std::vector<double>(ncz, 0.0)));
    p.resize(ncx, std::vector<std::vector<double>>(ncy, std::vector<double>(ncz, 0.0)));
    t.resize(ncx, std::vector<std::vector<double>>(ncy, std::vector<double>(ncz, 0.0)));

    // u0, v0, w0, p0, t0
    u0.resize(ncx, std::vector<std::vector<double>>(ncy, std::vector<double>(ncz, 0.0)));
    v0.resize(ncx, std::vector<std::vector<double>>(ncy, std::vector<double>(ncz, 0.0)));
    w0.resize(ncx, std::vector<std::vector<double>>(ncy, std::vector<double>(ncz, 0.0)));
    p0.resize(ncx, std::vector<std::vector<double>>(ncy, std::vector<double>(ncz, 0.0)));
    t0.resize(ncx, std::vector<std::vector<double>>(ncy, std::vector<double>(ncz, 0.0)));

    // uf, vf, wf
    initial_uf = 0.0;
    uf.resize(nx, std::vector<std::vector<double>>(ncy, std::vector<double>(ncz, 0.0)));
    initial_vf = 0.0;
    vf.resize(ncx, std::vector<std::vector<double>>(ny, std::vector<double>(ncz, 0.0)));
    initial_wf = 0.0;
    wf.resize(ncx, std::vector<std::vector<double>>(ncy, std::vector<double>(nz, 0.0)));

    // 压力修正 pp
    pp.resize(ncx, std::vector<std::vector<double>>(ncy, std::vector<double>(ncz, 0.0)));

    // 残差 ru, rv, rw, rp, rt
    ru.resize(nx, std::vector<std::vector<double>>(ncy, std::vector<double>(ncz, 0.0)));
    rv.resize(ncx, std::vector<std::vector<double>>(ny, std::vector<double>(ncz, 0.0)));
    rw.resize(ncx, std::vector<std::vector<double>>(ncy, std::vector<double>(nz, 0.0)));
    rp.resize(ncx, std::vector<std::vector<double>>(ncy, std::vector<double>(ncz, 0.0)));
    rt.resize(ncx, std::vector<std::vector<double>>(ncy, std::vector<double>(ncz, 0.0)));
}

// 初始化

void StructuredMesh::SetInitialU(double U) {
    for (int i = 0; i < ncx; ++i)
        for (int j = 0; j < ncy; ++j)
            for (int k = 0; k < ncz; ++k)
                u[i][j][k] = U;
}


void StructuredMesh::SetInitialV(double V) {
    for (int i = 0; i < ncx; ++i)
        for (int j = 0; j < ncy; ++j)
            for (int k = 0; k < ncz; ++k)
                v[i][j][k] = V;
}


void StructuredMesh::SetInitialW(double W) {
    for (int i = 0; i < ncx; ++i)
        for (int j = 0; j < ncy; ++j)
            for (int k = 0; k < ncz; ++k)
                w[i][j][k] = W;
}


void StructuredMesh::SetInitialP(double P) {
    for (int i = 0; i < ncx; ++i)
        for (int j = 0; j < ncy; ++j)
            for (int k = 0; k < ncz; ++k)
                p[i][j][k] = P;
}


void StructuredMesh::SetInitialT(double T) {
    for (int i = 0; i < ncx; ++i)
        for (int j = 0; j < ncy; ++j)
            for (int k = 0; k < ncz; ++k)
                t[i][j][k] = T;
}


void StructuredMesh::SetInitialUF(double uf_value) {
    initial_uf = uf_value;
    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ncy; ++j)
            for (int k = 0; k < ncz; ++k)
                uf[i][j][k] = uf_value;
}


void StructuredMesh::SetInitialVF(double vf_value) {
    initial_vf = vf_value;
    for (int i = 0; i < ncx; ++i)
        for (int j = 0; j < ny; ++j)
            for (int k = 0; k < ncz; ++k)
                vf[i][j][k] = vf_value;
}


void StructuredMesh::SetInitialWF(double wf_value) {
    initial_wf = wf_value;
    for (int i = 0; i < ncx; ++i)
        for (int j = 0; j < ncy; ++j)
            for (int k = 0; k < nz; ++k)
                wf[i][j][k] = wf_value;
}


double StructuredMesh::GetInitialUF() const {
    return initial_uf;
}


double StructuredMesh::GetInitialVF() const {
    return initial_vf;
}


double StructuredMesh::GetInitialWF() const {
    return initial_wf;
}


void StructuredMesh::CreateCoeffMeshData() {

    if (dim == 2) {
        ncoef = 6; // aP, aE, aW, aN, aS, bsrc
    } else {
        ncoef = 8; // aP, aE, aW, aN, aS, aT, aB, bsrc
    }
    ncoef_p = ncoef; // 用于分离求解器的压力系数矩阵

    // 设置系数索引
    id_aP = 0;
    id_aE = 1;
    id_aW = 2;
    id_aN = 3;
    id_aS = 4;
    if (dim == 3) {
        id_aT = 5;
        id_aB = 6;
    }

    // 源项系数
    id_bsrc = ncoef - 1;

    // 初始化各方向的系数矩阵 cu, cv, cw, cp, ct
    cu.resize(ncx, std::vector<std::vector<std::vector<double>>>(ncy,
        std::vector<std::vector<double>>(ncz, std::vector<double>(ncoef, 0.0))));
    cv.resize(ncx, std::vector<std::vector<std::vector<double>>>(ncy,
        std::vector<std::vector<double>>(ncz, std::vector<double>(ncoef, 0.0))));
    cw.resize(ncx, std::vector<std::vector<std::vector<double>>>(ncy,
        std::vector<std::vector<double>>(ncz, std::vector<double>(ncoef, 0.0))));
    cp.resize(ncx, std::vector<std::vector<std::vector<double>>>(ncy,
        std::vector<std::vector<double>>(ncz, std::vector<double>(ncoef_p, 0.0))));
    ct.resize(ncx, std::vector<std::vector<std::vector<double>>>(ncy,
        std::vector<std::vector<double>>(ncz, std::vector<double>(ncoef, 0.0))));
}


void StructuredMesh::CreateSimulationData() {

    grav.resize(dim, 0.0);

    // 标识符
    ieqn = 1;
    eqn_conduction = 0;
    eqn_flow = 1;
    eqn_conduction_flow = 2;

    // 仿真参数初始化
    nsteps = 1;
    dt = 1.0;
    stop_sim = false;

    conv_scheme = 0;

    imethod = 0;
    method_simple = 0;
    method_simplec = 1;

    isolver = 0;
    solver_pj = 0;      // 投影法
    solver_gs = 1;      // 高斯法

    itest = 0;
    test_lid = 0;
    test_channel = 1;
    test_channel_temp = 2;
}


void StructuredMesh::CreateSolvingMethodData() {
    
    // 最大迭代次数设置
    niter_u = 100;
    niter_v = 100;
    niter_w = 100;
    niter_p = 100;
    niter_t = 10;

    // 松弛因子
    relax_u = 0.0;
    relax_v = 0.0;
    relax_w = 0.0;
    relax_p = 0.0;
    relax_t = 0.75;

    // 残差阈值
    res_u = 1e-2;
    res_v = 1e-2;
    res_w = 1e-2;
    res_p = 1e-2;
    res_t = 1e-2;

    total_linsol_iters = 0;

     // 容差
    mass_tol = 1e-6;
    temp_tol = 1e-6;

    // 初始化
    l2_curr = 0.0;
    l2_max = -1e20;
    l2_max_u = -1e20;
    l2_max_v = -1e20;
    l2_max_w = -1e20;
    l2_max_p = -1e20;
    l2_max_pp = -1e20;
    l2_max_t = -1e20;
    l2_u = 0.0;
    l2_v = 0.0;
    l2_w = 0.0;
    l2_p = 0.0;
    l2_pp = 0.0;
    l2_t = 0.0;
}


void StructuredMesh::Set_ieqn(int ieqn) {
    this->ieqn = ieqn;
}


void StructuredMesh::Set_nsteps(int nsteps) {
    this->nsteps = nsteps;
}


void StructuredMesh::Set_dt(double dt) {
    this->dt = dt;
}


void StructuredMesh::Set_conv_scheme(int conv_scheme) {
    this->conv_scheme = conv_scheme;
}


void StructuredMesh::Set_temp_solver_param(int niter_t, double relax_t, double res_t, double temp_tol) {
    this->niter_t = niter_t;
    this->relax_t = relax_t;
    this->res_t = res_t;
    this->temp_tol = temp_tol;
}


void StructuredMesh::Set_xvel_solver_param(int niter_u, double relax_u, double res_u) {
    this->niter_u = niter_u;
    this->relax_u = relax_u;
    this->res_u = res_u;
}


void StructuredMesh::Set_yvel_solver_param(int niter_v, double relax_v, double res_v) {
    this->niter_v = niter_v;
    this->relax_v = relax_v;
    this->res_v = res_v;
}


void StructuredMesh::Set_zvel_solver_param(int niter_w, double relax_w, double res_w) {
    this->niter_w = niter_w;
    this->relax_w = relax_w;
    this->res_w = res_w;
}


void StructuredMesh::Set_pre_solver_param(int niter_p, double relax_p, double res_p) {
    this->niter_p = niter_p;
    this->relax_p = relax_p;
    this->res_p = res_p;
}


void StructuredMesh::Set_mass_tol(double mass_tol) {
    this->mass_tol = mass_tol;
}



Fluid::Fluid(int ncx, int ncy, int ncz)
    : ncx(ncx), ncy(ncy), ncz(ncz),con(0.0), spht(0.0), heat_src(0.0){

    // 初始化密度和粘度场，初始值为 1.0
    dens.resize(ncx, std::vector<std::vector<double>>(ncy, std::vector<double>(ncz, 1.0)));
    mu.resize(ncx, std::vector<std::vector<double>>(ncy, std::vector<double>(ncz, 1.0)));
}

// 设置初始密度和粘度
void Fluid::SetInitialDenMu(double dens_val, double mu_val) {
    for (int i = 0; i < ncx; ++i) {
        for (int j = 0; j < ncy; ++j) {
            for (int k = 0; k < ncz; ++k) {
                dens[i][j][k] = dens_val;
                mu[i][j][k] = mu_val;
            }
        }
    }
}

// 设置导热系数和比热容
void Fluid::SetConSpht(double con_val, double spht_val) {
    con = con_val;
    spht = spht_val;
}


BoundaryCondition::BoundaryCondition()
    : type(0) {}

// 速度边界
BoundaryConditionVel::BoundaryConditionVel()
    : u_type(0), v_type(0), w_type(0),
      u(0.0), v(0.0), w(0.0),
      u_flux(0.0), v_flux(0.0), w_flux(0.0) {}

// 温度边界
BoundaryConditionTemp::BoundaryConditionTemp()
    : temp_type(0), t(0.0), heat_flux(0.0) {}


// 流体边界构造
FluidBoundary::FluidBoundary(int dim)
    : dim(dim),
    // 各方向
    fid_e(0), fid_w(1), fid_n(2), fid_s(3), fid_t(4), fid_b(5),
    // 边界条件类型
    bc_none(0), bc_wall(1), bc_inlet(2), bc_outlet(3),
    // 边界条件编号
    bcid_none(0), bcid_xmin(1), bcid_xmax(2), bcid_ymin(3), bcid_ymax(4), bcid_zmin(5), bcid_zmax(6),
    //常量和通量
    num_bcs(7), bc_constant(0), bc_flux(1),
    // 温度
    temp_east(0.0), temp_west(0.0), temp_north(0.0), temp_south(0.0), temp_top(0.0), temp_bottom(0.0),
    // 压力
    p_outlet(0.0), pp_outlet(0.0){

    // 初始化边界条件向量
    bcs.resize(num_bcs);
    bcs_vel.resize(num_bcs);
    bcs_temp.resize(num_bcs);
}


//创建单元面的边界条件
void FluidBoundary::CreateBoundaryOfCellFaces(const StructuredMesh& caseData) {

    int ncx = caseData.GetNCX();
    int ncy = caseData.GetNCY();
    int ncz = caseData.GetNCZ();

    const std::vector<double>& x = caseData.GetX();
    const std::vector<double>& y = caseData.GetY();
    const std::vector<double>& z = caseData.GetZ();

    double xmin = caseData.GetXMin();
    double xmax = caseData.GetXMax();
    double ymin = caseData.GetYMin();
    double ymax = caseData.GetYMax();
    double zmin = caseData.GetZMin();
    double zmax = caseData.GetZMax();

    // 初始化边界条件 ID 的三维数组
    bcid.resize(ncx, std::vector<std::vector<std::vector<int>>>(ncy,
        std::vector<std::vector<int>>(ncz, std::vector<int>(2 * dim, bcid_none))));

    double eps = 1e-12;

    // 遍历每个网格单元，设置边界条件
    for (int k = 0; k < ncz; ++k) {
        for (int j = 0; j < ncy; ++j) {
            for (int i = 0; i < ncx; ++i) {
                // X-direction boundaries
                double x0 = x[i];
                if (std::abs(x0 - xmin) < eps) {
                    bcid[i][j][k][fid_w] = bcid_xmin;
                }
                x0 = x[i + 1];
                if (std::abs(x0 - xmax) < eps) {
                    bcid[i][j][k][fid_e] = bcid_xmax;
                }

                // Y-direction boundaries
                double y0 = y[j];
                if (std::abs(y0 - ymin) < eps) {
                    bcid[i][j][k][fid_s] = bcid_ymin;
                }
                y0 = y[j + 1];
                if (std::abs(y0 - ymax) < eps) {
                    bcid[i][j][k][fid_n] = bcid_ymax;
                }

                if (dim == 3) {
                    // Z-direction boundaries
                    double z0 = z[k];
                    if (std::abs(z0 - zmin) < eps) {
                        bcid[i][j][k][fid_b] = bcid_zmin;
                    }
                    z0 = z[k + 1];
                    if (std::abs(z0 - zmax) < eps) {
                        bcid[i][j][k][fid_t] = bcid_zmax;
                    }
                }
            }
        }
    }
}

// 创建边界条件
void FluidBoundary::CreateBoundaryData(int dim,
                                       const std::string& input_bc_xmin, const std::string& input_bc_xmax,
                                       const std::string& input_bc_ymin, const std::string& input_bc_ymax,
                                       const std::string& input_bc_zmin, const std::string& input_bc_zmax) {
    
    int id = bcid_none;

    // 默认为空
    bcs[id].type = bc_none;

    // 以下设置各个边界
    id = bcid_xmin;
    if (input_bc_xmin == "inlet") {
        bcs[id].type = bc_inlet;
    } else if (input_bc_xmin == "outlet") {
        bcs[id].type = bc_outlet;
    } else if (input_bc_xmin == "wall") {
        bcs[id].type = bc_wall;
    }


    id = bcid_xmax;
    if (input_bc_xmax == "inlet") {
        bcs[id].type = bc_inlet;
    } else if (input_bc_xmax == "outlet") {
        bcs[id].type = bc_outlet;
    } else if (input_bc_xmax == "wall") {
        bcs[id].type = bc_wall;
    }


    id = bcid_ymin;
    if (input_bc_ymin == "inlet") {
        bcs[id].type = bc_inlet;
    } else if (input_bc_ymin == "outlet") {
        bcs[id].type = bc_outlet;
    } else if (input_bc_ymin == "wall") {
        bcs[id].type = bc_wall;
    }


    id = bcid_ymax;
    if (input_bc_ymax == "inlet") {
        bcs[id].type = bc_inlet;
    } else if (input_bc_ymax == "outlet") {
        bcs[id].type = bc_outlet;
    } else if (input_bc_ymax == "wall") {
        bcs[id].type = bc_wall;
    }

    if (dim == 3) {

        id = bcid_zmin;
        if (input_bc_zmin == "inlet") {
            bcs[id].type = bc_inlet;
        } else if (input_bc_zmin == "outlet") {
            bcs[id].type = bc_outlet;
        } else if (input_bc_zmin == "wall") {
            bcs[id].type = bc_wall;
        }

        id = bcid_zmax;
        if (input_bc_zmax == "inlet") {
            bcs[id].type = bc_inlet;
        } else if (input_bc_zmax == "outlet") {
            bcs[id].type = bc_outlet;
        } else if (input_bc_zmax == "wall") {
            bcs[id].type = bc_wall;
        }
    }
}

// 速度边界条件
void FluidBoundary::CreateBoundaryDataVel(int dim,
    const std::string& input_bc_xmin, const std::vector<std::string>& input_vel_type_xmin, const std::vector<double>& input_vel_xmin,
    const std::string& input_bc_xmax, const std::vector<std::string>& input_vel_type_xmax, const std::vector<double>& input_vel_xmax,
    const std::string& input_bc_ymin, const std::vector<std::string>& input_vel_type_ymin, const std::vector<double>& input_vel_ymin,
    const std::string& input_bc_ymax, const std::vector<std::string>& input_vel_type_ymax, const std::vector<double>& input_vel_ymax,
    const std::string& input_bc_zmin, const std::vector<std::string>& input_vel_type_zmin,
    const std::vector<double>& input_vel_zmin,
    const std::string& input_bc_zmax, const std::vector<std::string>& input_vel_type_zmax,
    const std::vector<double>& input_vel_zmax) {

    int id = bcid_none;
    bcs_vel[id].u = 0.0;
    bcs_vel[id].v = 0.0;
    bcs_vel[id].u_flux = 0.0;
    bcs_vel[id].v_flux = 0.0;
    bcs_vel[id].u_type = bc_none;
    bcs_vel[id].v_type = bc_none;
    if (dim == 3) {
        bcs_vel[id].w = 0.0;
        bcs_vel[id].w_flux = 0.0;
        bcs_vel[id].w_type = bc_none;
    }

    auto set_vel_bc = [this](int id, const std::vector<std::string>& vel_types, const std::vector<double>& vel_values, int dim) {

        if (vel_types.size() > 0) {
            if (vel_types[0] == "constant") {
                bcs_vel[id].u_type = bc_constant;
                bcs_vel[id].u = vel_values[0];
            } else if (vel_types[0] == "flux") {
                bcs_vel[id].u_type = bc_flux;
                bcs_vel[id].u_flux = vel_values[0];
            }
        }

        if (vel_types.size() > 1) {
            if (vel_types[1] == "constant") {
                bcs_vel[id].v_type = bc_constant;
                bcs_vel[id].v = vel_values[1];
            } else if (vel_types[1] == "flux") {
                bcs_vel[id].v_type = bc_flux;
                bcs_vel[id].v_flux = vel_values[1];
            }
        }

        if (dim == 3 && vel_types.size() > 2) {
            if (vel_types[2] == "constant") {
                bcs_vel[id].w_type = bc_constant;
                bcs_vel[id].w = vel_values[2];
            } else if (vel_types[2] == "flux") {
                bcs_vel[id].w_type = bc_flux;
                bcs_vel[id].w_flux = vel_values[2];
            }
        }
    };

    id = bcid_xmin;
    set_vel_bc(id, input_vel_type_xmin, input_vel_xmin, dim);

    id = bcid_xmax;
    set_vel_bc(id, input_vel_type_xmax, input_vel_xmax, dim);

    id = bcid_ymin;
    set_vel_bc(id, input_vel_type_ymin, input_vel_ymin, dim);

    id = bcid_ymax;
    set_vel_bc(id, input_vel_type_ymax, input_vel_ymax, dim);

    if (dim == 3) {
        id = bcid_zmin;
        set_vel_bc(id, input_vel_type_zmin, input_vel_zmin, dim);

        id = bcid_zmax;
        set_vel_bc(id, input_vel_type_zmax, input_vel_zmax, dim);
    }
}

// 温度边界条件
void FluidBoundary::CreateBoundaryDataTemp(int dim,
                                           const std::string& input_bc_xmin, const std::string& input_temp_type_xmin, double input_heat_flux_xmin,
                                           const std::string& input_bc_xmax, const std::string& input_temp_type_xmax, double input_heat_flux_xmax,
                                           const std::string& input_bc_ymin, const std::string& input_temp_type_ymin, double input_heat_flux_ymin,
                                           const std::string& input_bc_ymax, const std::string& input_temp_type_ymax, double input_heat_flux_ymax) {
    int id = bcid_none;
    bcs_temp[id].t = 0.0;
    bcs_temp[id].heat_flux = 0.0;
    bcs_temp[id].temp_type = bc_none;

    auto set_temp_bc = [this](int id, const std::string& temp_type, double value) {
        if (temp_type == "constant") {
            bcs_temp[id].temp_type = bc_constant;
            bcs_temp[id].t = value;
        } else if (temp_type == "heat_flux") {
            bcs_temp[id].temp_type = bc_flux;
            bcs_temp[id].heat_flux = value;
        }
    };

    id = bcid_xmin;
    set_temp_bc(id, input_temp_type_xmin, input_heat_flux_xmin);

    id = bcid_xmax;
    set_temp_bc(id, input_temp_type_xmax, input_heat_flux_xmax);

    id = bcid_ymin;
    set_temp_bc(id, input_temp_type_ymin, input_heat_flux_ymin);

    id = bcid_ymax;
    set_temp_bc(id, input_temp_type_ymax, input_heat_flux_ymax);
}
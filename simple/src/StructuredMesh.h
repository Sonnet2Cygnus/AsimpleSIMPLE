#ifndef STRUCTUREDMESH_H
#define STRUCTUREDMESH_H

#include <vector>
#include <string>
#include <iostream>
#include <cmath>


// double 精度，与python 中精度 fp 相对应
using fp = double;


// 前向声明，在编译前提前告诉编译器类的种类与位置，适合多种类复合使用的时候
class Fluid;
class BoundaryCondition;
class BoundaryConditionVel;
class BoundaryConditionTemp;
class FluidBoundary;


class StructuredMesh {
public:
    
    // 构造函数
    StructuredMesh();

    // 网格生成
    void CreateMesh(int dim, int ncx, int ncy, int ncz = 1);                        // 指定维度和网格单元数
    void CreateCoordinates(double xmin, double xmax, double ymin, double ymax,      // 定义网格坐标范围
                           double zmin = 0.0, double zmax = 1.0);         

    // 场数据，系数矩阵，仿真数据和算法数据          
    void CreateFieldMeshData();                                                    
    void CreateCoeffMeshData();
    void CreateSimulationData();
    void CreateSolvingMethodData();

    // 定义初始条件
    void SetInitialU(double U);
    void SetInitialV(double V);
    void SetInitialW(double W);
    void SetInitialP(double P);
    void SetInitialT(double T);
    // 面速度
    void SetInitialUF(double uf_value);
    void SetInitialVF(double vf_value);
    void SetInitialWF(double wf_value);

    // 获取初始面速度
    double GetInitialUF() const;
    double GetInitialVF() const;
    double GetInitialWF() const;

    // 设置仿真参数
    void Set_ieqn(int ieqn);
    void Set_nsteps(int nsteps);
    void Set_dt(double dt);
    void Set_conv_scheme(int conv_scheme);

    // 设置求解器参数
    void Set_temp_solver_param(int niter_t, double relax_t, double res_t, double temp_tol);
    void Set_xvel_solver_param(int niter_u, double relax_u, double res_u);
    void Set_yvel_solver_param(int niter_v, double relax_v, double res_v);
    void Set_zvel_solver_param(int niter_w, double relax_w, double res_w);
    void Set_pre_solver_param(int niter_p, double relax_p, double res_p);
    void Set_mass_tol(double mass_tol);

    // 获取网格参数的函数
    int GetNCX() const { return ncx; }
    int GetNCY() const { return ncy; }
    int GetNCZ() const { return ncz; }
    const std::vector<double>& GetX() const { return x; }
    const std::vector<double>& GetY() const { return y; }
    const std::vector<double>& GetZ() const { return z; }
    double GetXMin() const { return xmin; }
    double GetXMax() const { return xmax; }
    double GetYMin() const { return ymin; }
    double GetYMax() const { return ymax; }
    double GetZMin() const { return zmin; }
    double GetZMax() const { return zmax; }

    // 网格变量
    int dim;
    int ncx, ncy, ncz;
    int nx, ny, nz;
    double dx, dy, dz;
    std::vector<double> x, y, z;
    std::vector<double> xc, yc, zc;

    // 场变量
    std::vector<std::vector<std::vector<double>>> u, v, w, p, t;
    std::vector<std::vector<std::vector<double>>> u0, v0, w0, p0, t0;

    // 面变量
    double initial_uf, initial_vf, initial_wf;
    std::vector<std::vector<std::vector<double>>> uf, vf, wf;
    std::vector<std::vector<std::vector<double>>> pp;

    // 残差
    std::vector<std::vector<std::vector<double>>> ru, rv, rw, rp, rt;

    // 系数矩阵
    int ncoef, ncoef_p;
    int id_aP, id_aE, id_aW, id_aN, id_aS, id_aT, id_aB, id_bsrc;
    std::vector<std::vector<std::vector<std::vector<double>>>> cu, cv, cw, cp, ct;

    // 仿真参数
    std::vector<double> grav;
    int ieqn, eqn_conduction, eqn_flow, eqn_conduction_flow;
    int nsteps;
    double dt;
    bool stop_sim;
    int conv_scheme;
    int imethod, method_simple, method_simplec;                     // 求解方法选择
    int isolver, solver_pj, solver_gs;                              // 求解器类型
    int itest, test_lid, test_channel, test_channel_temp;           // 测试场景标识

    // 求解器参数
    int niter_u, niter_v, niter_w, niter_p, niter_t;                // 最大迭代次数
    double relax_u, relax_v, relax_w, relax_p, relax_t;             // 松弛因子
    double res_u, res_v, res_w, res_p, res_t;
    int total_linsol_iters;                                         // 线性求解器迭代次数
    double mass_tol, temp_tol;
    double l2_curr, l2_max, l2_max_u, l2_max_v, l2_max_w, l2_max_p, l2_max_pp, l2_max_t;
    double l2_u, l2_v, l2_w, l2_p, l2_pp, l2_t;

   // 边界条件标识符
    int fid_e; 
    int fid_w; 
    int fid_n; 
    int fid_s; 
    int fid_t; 
    int fid_b; 

    // 边界条件类型
    int bc_none;
    int bc_wall;
    int bc_inlet;
    int bc_outlet;

    // 边界条件ID
    int bcid_none;
    int bcid_xmin;
    int bcid_xmax;
    int bcid_ymin;
    int bcid_ymax;
    int bcid_zmin;
    int bcid_zmax;

    // 边界条件数量
    int num_bcs;

    // 边界条件 常数。通量
    int bc_constant;
    int bc_flux;

    // 温度边界条件
    double temp_east, temp_west, temp_north, temp_south, temp_top, temp_bottom;

    // 压力边界条件
    double p_outlet, pp_outlet;

    // 密度和粘度场
    std::vector<std::vector<std::vector<double>>> dens, mu;

    // 边界条件列表
    std::vector<BoundaryCondition> bcs;
    std::vector<BoundaryConditionVel> bcs_vel;
    std::vector<BoundaryConditionTemp> bcs_temp;

    // 每个单元面的边界条件ID
    std::vector<std::vector<std::vector<std::vector<int>>>> bcid;

private:
    // 网格边界坐标
    double xmin, xmax, ymin, ymax, zmin, zmax;
};


class Fluid {
public:
    // 构造函数
    Fluid(int ncx, int ncy, int ncz);

    // 密度；粘度；导热系数和比热
    void SetInitialDenMu(double dens_val, double mu_val);
    void SetConSpht(double con_val, double spht_val);

    //网格尺寸
    int ncx, ncy, ncz;

    //密度和粘度场
    std::vector<std::vector<std::vector<double>>> dens, mu;

    // 导热系数和比热
    double con, spht;

    // 热源
    double heat_src;
};


class BoundaryCondition {
public:
    BoundaryCondition();
    int type;
};


// 速度边界
class BoundaryConditionVel {
public:
    BoundaryConditionVel();
    int u_type, v_type, w_type;
    double u, v, w;
    double u_flux, v_flux, w_flux;
};


// 温度边界
class BoundaryConditionTemp {
public:
    BoundaryConditionTemp();
    int temp_type;
    double t;
    double heat_flux;
};


// 流体边界
class FluidBoundary {
public:
     // 构造函数
    FluidBoundary(int dim);

    // 根据单元面创建边界
    void CreateBoundaryOfCellFaces(const StructuredMesh& caseData);
    void CreateBoundaryData(int dim,
                            const std::string& input_bc_xmin, const std::string& input_bc_xmax,
                            const std::string& input_bc_ymin, const std::string& input_bc_ymax,
                            const std::string& input_bc_zmin = "wall", const std::string& input_bc_zmax = "wall");
    
    // 创建速度边界数据
    void CreateBoundaryDataVel(int dim,
                               const std::string& input_bc_xmin, const std::vector<std::string>& input_vel_type_xmin, const std::vector<double>& input_vel_xmin,
                               const std::string& input_bc_xmax, const std::vector<std::string>& input_vel_type_xmax, const std::vector<double>& input_vel_xmax,
                               const std::string& input_bc_ymin, const std::vector<std::string>& input_vel_type_ymin, const std::vector<double>& input_vel_ymin,
                               const std::string& input_bc_ymax, const std::vector<std::string>& input_vel_type_ymax, const std::vector<double>& input_vel_ymax,
                               const std::string& input_bc_zmin = "wall", const std::vector<std::string>& input_vel_type_zmin = {},
                               const std::vector<double>& input_vel_zmin = {},
                               const std::string& input_bc_zmax = "wall", const std::vector<std::string>& input_vel_type_zmax = {},
                               const std::vector<double>& input_vel_zmax = {});

    // 创建温度边界数据
    void CreateBoundaryDataTemp(int dim,
                                const std::string& input_bc_xmin, const std::string& input_temp_type_xmin, double input_heat_flux_xmin,
                                const std::string& input_bc_xmax, const std::string& input_temp_type_xmax, double input_heat_flux_xmax,
                                const std::string& input_bc_ymin, const std::string& input_temp_type_ymin, double input_heat_flux_ymin,
                                const std::string& input_bc_ymax, const std::string& input_temp_type_ymax, double input_heat_flux_ymax);

    // 变量
    int fid_e, fid_w, fid_n, fid_s, fid_t, fid_b;                                           //各个面标识符
    int bc_none, bc_wall, bc_inlet, bc_outlet;
    int bcid_none, bcid_xmin, bcid_xmax, bcid_ymin, bcid_ymax, bcid_zmin, bcid_zmax;        // 边界条件编号
    int num_bcs, bc_constant, bc_flux;                                                      
    double temp_east, temp_west, temp_north, temp_south, temp_top, temp_bottom;
    double p_outlet, pp_outlet;

    // 边界列表，速度和温度
    std::vector<BoundaryCondition> bcs;                                                    
    std::vector<BoundaryConditionVel> bcs_vel;
    std::vector<BoundaryConditionTemp> bcs_temp;

     // 每个单元面的边界条件编号
    std::vector<std::vector<std::vector<std::vector<int>>>> bcid;

private:
    int dim;
};

#endif 

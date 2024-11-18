#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include "StructuredMesh.h"
#include "collocated_segregated.h"

int main() {

    // 开始计时
    auto clock_begin = std::chrono::high_resolution_clock::now();

    // 设置维度和网格信息
    int dim = 2;
    int ncx = 32, ncy = 64, ncz = 1;

    double xmin = -0.0625, xmax = 0.0625;
    double ymin = -0.0125, ymax = 0.0125;
    double zmin = 0.0, zmax = 1.0;

    // 创建网格对象
    StructuredMesh caseData;
    caseData.CreateMesh(dim, ncx, ncy, ncz);
    caseData.CreateCoordinates(xmin, xmax, ymin, ymax, zmin, zmax);
    caseData.CreateFieldMeshData();

    // 设置初始温度场
    double tref = 288.0;
    caseData.SetInitialT(tref);
    caseData.CreateCoeffMeshData();

    // 设置模拟参数
    caseData.CreateSimulationData();
    int conv_scheme = 0;
    caseData.Set_conv_scheme(conv_scheme);

    // 选择测试方程
    int ieqn = 2;
    caseData.Set_ieqn(ieqn);

    // 设置迭代次数
    int nsteps = 800;
    caseData.Set_nsteps(nsteps);

    double dt = 1.e16;
    caseData.Set_dt(dt);
    caseData.CreateSolvingMethodData();


    // 温度求解器参数
    int niter_t = 200;
    double relax_t = 0.75;
    double res_t = 1.e-4;
    double temp_tol = 1.e-6;
    caseData.Set_temp_solver_param(niter_t, relax_t, res_t, temp_tol);


    // 设置速度和压力求解器参数
    caseData.Set_xvel_solver_param(5, 0.75, 1.e-4);
    caseData.Set_yvel_solver_param(5, 0.75, 1.e-4);
    caseData.Set_pre_solver_param(200, 0.25, 1.e-6);
    caseData.Set_mass_tol(1.e-6);


    // 创建流体并设置属性
    Fluid fluid(ncx, ncy, ncz);
    double dens = 1.0, mu = 0.01;
    fluid.SetInitialDenMu(dens, mu);

    double con = 1.0, spht = 1.0;
    fluid.SetConSpht(con, spht);


    // 创建边界条件
    FluidBoundary fluidBoundary(dim);
    fluidBoundary.CreateBoundaryOfCellFaces(caseData);
    fluidBoundary.CreateBoundaryData(dim, "inlet", "outlet", "wall", "wall", "wall", "wall");
    fluidBoundary.CreateBoundaryDataVel(dim,
        "inlet", {"constant", "constant", "flux"}, {1.0, 0.0, 0.0},
        "outlet", {"None", "None", "None"}, {0.0, 0.0, 0.0},
        "wall", {"constant", "constant", "flux"}, {0.0, 0.0, 0.0},
        "wall", {"constant", "constant", "flux"}, {0.0, 0.0, 0.0},
        "wall", {"constant", "constant", "flux"}, {0.0, 0.0, 0.0},
        "wall", {"constant", "constant", "flux"}, {0.0, 0.0, 0.0});
    fluidBoundary.CreateBoundaryDataTemp(dim,
        "inlet", "constant", 300.0,
        "outlet", "heat_flux", 0.0,
        "wall", "constant", 320.0,
        "wall", "heat_flux", 0.0);

    // 主要迭代流程
    for (int it = 1; it <= nsteps; ++it) {
        if (it % 1 == 0 || it == 1 || it == nsteps) {
            std::cout << "\n----------------------------\n";
            std::cout << "Begin iter = " << it << "\n";
            std::cout << "----------------------------\n";
        }

        caseData.stop_sim = collocated_segregated(it, caseData, fluid, fluidBoundary);

        if (caseData.stop_sim) break;

        if (it % 1 == 0 || it == nsteps || caseData.stop_sim) {
            std::cout << "it: " << it << " walltime: "
                      << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - clock_begin).count()
                      << " l2_u/l2_max_u: " << caseData.l2_u / caseData.l2_max_u
                      << " l2_v/l2_max_v: " << caseData.l2_v / caseData.l2_max_v
                      << " l2_p/l2_max_p: " << caseData.l2_p / caseData.l2_max_p
                      << " l2_pp/l2_max_pp: " << caseData.l2_pp / caseData.l2_max_pp
                      << "\n";
        }
    }

    // 输出总耗时
    auto elapsed_time = std::chrono::high_resolution_clock::now() - clock_begin;
    std::cout << "Elapsed time: " << std::chrono::duration<double>(elapsed_time).count() << " seconds\n";

    return 0;
}

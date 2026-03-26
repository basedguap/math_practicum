#include <gmsh.h>
#include <cmath>
#include <vector>

int main()
{
    gmsh::initialize();
    gmsh::model::add("hollow_torus");

    const double pi = 3.141592653589793;
    const double R = 5.0;
    const double r_outer = 1.5;
    const double r_inner = 1.0;
    const double lc = (r_outer - r_inner) / 4.0;
    const int num_toroidal = 16;
    const int num_poloidal = 8;

    std::vector<std::vector<int>> pt_outer(num_toroidal, std::vector<int>(num_poloidal));
    std::vector<std::vector<int>> pt_inner(num_toroidal, std::vector<int>(num_poloidal));
    std::vector<int> center_pol(num_toroidal);
    std::vector<int> center_tor_outer(num_poloidal);
    std::vector<int> center_tor_inner(num_poloidal);

    for (int i = 0; i < num_toroidal; i++) {
        double u = 2.0 * pi * i / num_toroidal;
        center_pol[i] = gmsh::model::geo::addPoint(R * cos(u), R * sin(u), 0, lc);
        for (int j = 0; j < num_poloidal; j++) {
            double v = 2.0 * pi * j / num_poloidal;
            pt_outer[i][j] = gmsh::model::geo::addPoint((R + r_outer * cos(v)) * cos(u), (R + r_outer * cos(v)) * sin(u), r_outer * sin(v), lc);
            pt_inner[i][j] = gmsh::model::geo::addPoint((R + r_inner * cos(v)) * cos(u), (R + r_inner * cos(v)) * sin(u), r_inner * sin(v), lc);
        }
    }

    for (int j = 0; j < num_poloidal; j++) {
        double v = 2.0 * pi * j / num_poloidal;
        center_tor_outer[j] = gmsh::model::geo::addPoint(0, 0, r_outer * sin(v), lc);
        center_tor_inner[j] = gmsh::model::geo::addPoint(0, 0, r_inner * sin(v), lc);
    }

    std::vector<std::vector<int>> arc_pol_outer(num_toroidal, std::vector<int>(num_poloidal));
    std::vector<std::vector<int>> arc_pol_inner(num_toroidal, std::vector<int>(num_poloidal));
    std::vector<std::vector<int>> arc_tor_outer(num_toroidal, std::vector<int>(num_poloidal));
    std::vector<std::vector<int>> arc_tor_inner(num_toroidal, std::vector<int>(num_poloidal));

    for (int i = 0; i < num_toroidal; i++)
        for (int j = 0; j < num_poloidal; j++) {
            arc_pol_outer[i][j] = gmsh::model::geo::addCircleArc(pt_outer[i][j], center_pol[i], pt_outer[i][(j+1)%num_poloidal]);
            arc_pol_inner[i][j] = gmsh::model::geo::addCircleArc(pt_inner[i][j], center_pol[i], pt_inner[i][(j+1)%num_poloidal]);
            arc_tor_outer[i][j] = gmsh::model::geo::addCircleArc(pt_outer[i][j], center_tor_outer[j], pt_outer[(i+1)%num_toroidal][j]);
            arc_tor_inner[i][j] = gmsh::model::geo::addCircleArc(pt_inner[i][j], center_tor_inner[j], pt_inner[(i+1)%num_toroidal][j]);
        }

    std::vector<int> surf_outer, surf_inner;

    for (int i = 0; i < num_toroidal; i++)
        for (int j = 0; j < num_poloidal; j++) {
            std::vector<int> curves_outer(4);
            curves_outer[0] = arc_pol_outer[i][j];
            curves_outer[1] = arc_tor_outer[i][(j+1)%num_poloidal];
            curves_outer[2] = -arc_pol_outer[(i+1)%num_toroidal][j];
            curves_outer[3] = -arc_tor_outer[i][j];
            int loop_outer = gmsh::model::geo::addCurveLoop(curves_outer);

            std::vector<int> fill_outer(1);
            fill_outer[0] = loop_outer;
            surf_outer.push_back(gmsh::model::geo::addSurfaceFilling(fill_outer));

            std::vector<int> curves_inner(4);
            curves_inner[0] = arc_pol_inner[i][j];
            curves_inner[1] = arc_tor_inner[i][(j+1)%num_poloidal];
            curves_inner[2] = -arc_pol_inner[(i+1)%num_toroidal][j];
            curves_inner[3] = -arc_tor_inner[i][j];
            int loop_inner = gmsh::model::geo::addCurveLoop(curves_inner);

            std::vector<int> fill_inner(1);
            fill_inner[0] = loop_inner;
            surf_inner.push_back(gmsh::model::geo::addSurfaceFilling(fill_inner));
        }

    int shell_outer = gmsh::model::geo::addSurfaceLoop(surf_outer);
    int shell_inner = gmsh::model::geo::addSurfaceLoop(surf_inner);

    std::vector<int> shells(2);
    shells[0] = shell_outer;
    shells[1] = shell_inner;
    gmsh::model::geo::addVolume(shells);

    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate(3);
    gmsh::write("hollow_torus.msh");
    gmsh::fltk::run();
    gmsh::finalize();
    return 0;
}
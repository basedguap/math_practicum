#include <gmsh.h>
#include <cmath>
#include <vector>

int main() {
    gmsh::initialize();
    gmsh::model::add("pokemon");

    gmsh::merge("../char.stl");

    gmsh::model::mesh::classifySurfaces(M_PI * 40 / 180, true, false, M_PI);
    gmsh::model::mesh::createGeometry();

    std::vector<std::pair<int, int>> surfaces;
    gmsh::model::getEntities(surfaces, 2);

    std::vector<int> surface_tags;

    for (int i = 0; i < surfaces.size(); i++) {
        int current_tag = surfaces[i].second;
        surface_tags.push_back(current_tag);
        gmsh::model::setColor({{2, current_tag}}, 255, 0, 0);
    }

    gmsh::model::geo::synchronize();
    int surface_loop = gmsh::model::geo::addSurfaceLoop(surface_tags);
    gmsh::model::geo::addVolume({surface_loop});
    gmsh::model::geo::synchronize();

    gmsh::model::mesh::generate(3);

    gmsh::option::setColor("Mesh.Lines", 0, 0, 0);

    gmsh::write("Gengar.msh");

    gmsh::fltk::run();

    gmsh::finalize();

    return 0;
}
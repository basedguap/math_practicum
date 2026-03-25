#include <gmsh.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnsignedCharArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

int main () {
  gmsh::initialize();
  gmsh::model::add("anim");

  gmsh::merge("../char.stl");

  gmsh::model::mesh::removeDuplicateNodes();
  gmsh::model::mesh::classifySurfaces(40 * M_PI / 180, true, true, M_PI);
  gmsh::model::mesh::createGeometry();

  std::vector<std::pair<int, int>> surfaces;
  gmsh::model::getEntities(surfaces, 2);

  std::vector<int> surfaceTags;
  for (auto s : surfaces) surfaceTags.push_back(s.second);

  gmsh::model::geo::synchronize();
  int loop = gmsh::model::geo::addSurfaceLoop(surfaceTags);
  gmsh::model::geo::addVolume({loop});
  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate(3);

  std::vector<std::size_t> nodeTags;
  std::vector<double> nodeCoords, tmp;
  gmsh::model::mesh::getNodes(nodeTags, nodeCoords, tmp);

  std::vector<int> elementTypes;
  std::vector<std::vector<std::size_t>> elementTags, elementNodeTags;
  gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags, 3);

  int tetraIndex = -1;
  for (int i = 0; i < elementTypes.size(); i++) {
    if (elementTypes[i] == 4) { tetraIndex = i; break; }
  }

  int pointCount = nodeTags.size();

  std::size_t maxTag = 0;
  for (int i = 0; i < pointCount; i++) {
    if (nodeTags[i] > maxTag) maxTag = nodeTags[i];
  }

  std::vector<int> tagToIndex(maxTag + 1, -1);
  for (int i = 0; i < pointCount; i++) tagToIndex[nodeTags[i]] = i;

  double minX = 1e9, maxX = -1e9;
  double minY = 1e9, maxY = -1e9;
  double minZ = 1e9, maxZ = -1e9;

  for (int i = 0; i < pointCount; i++) {
    double x = nodeCoords[3 * i];
    double y = nodeCoords[3 * i + 1];
    double z = nodeCoords[3 * i + 2];

    if (x < minX) minX = x; if (x > maxX) maxX = x;
    if (y < minY) minY = y; if (y > maxY) maxY = y;
    if (z < minZ) minZ = z; if (z > maxZ) maxZ = z;
  }

  double centerX = (minX + maxX) / 2;
  double centerY = (minY + maxY) / 2;
  double centerZ = (minZ + maxZ) / 2;

  double sizeX = maxX - minX;
  double sizeY = maxY - minY;
  double sizeZ = maxZ - minZ;
  double maxSize = std::max(sizeX, std::max(sizeY, sizeZ));


  double tailTipX = 0, tailTipY = 0, tailTipZ = 0;
  double maxTailY = -1e9;
  for (int i = 0; i < pointCount; i++) {
    double x = nodeCoords[3 * i];
    double y = nodeCoords[3 * i + 1];
    double z = nodeCoords[3 * i + 2];
    if (z < minZ + 0.25 * sizeZ) continue;
    if (z > minZ + 0.85 * sizeZ) continue;
    if (y > maxTailY) {
      maxTailY = y;
      tailTipX = x; tailTipY = y; tailTipZ = z;
    }
  }


  double tailBaseX = centerX + 0.7 * (tailTipX - centerX);
  double tailBaseY = centerY + 0.7 * (tailTipY - centerY);
  double tailBaseZ = centerZ + 0.7 * (tailTipZ - centerZ);

  double tailDirX = tailTipX - tailBaseX;
  double tailDirY = tailTipY - tailBaseY;
  double tailDirZ = tailTipZ - tailBaseZ;
  double tdLen = sqrt(tailDirX*tailDirX + tailDirY*tailDirY + tailDirZ*tailDirZ);
  if (tdLen > 1e-6) { tailDirX /= tdLen; tailDirY /= tdLen; tailDirZ /= tdLen; }

  double tailRadius = tdLen * 1.8;

  int frames = 240;
  double totalTime = 3.0;
  double omega = 2 * M_PI * 1.5;
  double amplitude = 0.08 * sizeY;

  system("mkdir -p output");

  for (int f = 0; f < frames; f++) {
    double t = totalTime * f / (frames - 1);

    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(pointCount);

    auto velocity = vtkSmartPointer<vtkDoubleArray>::New();
    velocity->SetName("Velocity");
    velocity->SetNumberOfComponents(3);
    velocity->SetNumberOfTuples(pointCount);

    auto wave = vtkSmartPointer<vtkDoubleArray>::New();
    wave->SetName("Wave");
    wave->SetNumberOfTuples(pointCount);

    auto color = vtkSmartPointer<vtkUnsignedCharArray>::New();
    color->SetName("Color");
    color->SetNumberOfComponents(3);
    color->SetNumberOfTuples(pointCount);

    for (int i = 0; i < pointCount; i++) {
      double x = nodeCoords[3 * i];
      double y = nodeCoords[3 * i + 1];
      double z = nodeCoords[3 * i + 2];


      double shiftX = 0;
      double shiftY = 0;

      double relPosX = x - tailBaseX;
      double relPosY = y - tailBaseY;
      double relPosZ = z - tailBaseZ;

      double alongTail = relPosX * tailDirX + relPosY * tailDirY + relPosZ * tailDirZ;

      double perpX = relPosX - alongTail * tailDirX;
      double perpY = relPosY - alongTail * tailDirY;
      double perpZ = relPosZ - alongTail * tailDirZ;
      double pDist = sqrt(perpX*perpX + perpY*perpY + perpZ*perpZ);

      double weight = 0.0;

      if (alongTail > 0 && alongTail < tdLen * 2.0) {
        double skin = alongTail / tdLen;
        if (skin > 1.0) skin = 1.0;
        double blend = skin * skin * (3.0 - 2.0 * skin);
        double maxR = tailRadius * (0.3 + 0.7 * (1.0 - skin));
        double radFade = (pDist <= maxR) ? (1.0 - (pDist/maxR)*(pDist/maxR)) : 0.0;
        weight = blend * radFade;
        if (weight < 0) weight = 0;
        if (weight > 1) weight = 1;
      }

      if (weight > 0) {
        double angle = 0.5 * weight * sin(omega * t);
        double dx = x - tailBaseX;
        double dy = y - tailBaseY;
        double cosA = cos(angle);
        double sinA = sin(angle);
        double newX = tailBaseX + dx * cosA - dy * sinA;
        double newY = tailBaseY + dx * sinA + dy * cosA;
        shiftX = weight * (newX - x);
        shiftY = weight * (newY - y);
      }

      points->SetPoint(i, x + shiftX, y + shiftY, z);

      double vx = 0;
      double vy = 2 * amplitude * weight * omega * cos(omega * t); // ИЗМЕНЕНО: weight вместо линейного порога
      velocity->SetTuple3(i, vx, vy, 0);

      double radius = sqrt((x - centerX) * (x - centerX) +
                           (y - centerY) * (y - centerY) +
                           (z - centerZ) * (z - centerZ));

      double scalar = sin(8 * radius / (maxSize + 1e-12) - 2 * omega * t) +
                      0.4 * cos(6 * z / (sizeZ + 1e-12) + omega * t);

      wave->SetValue(i, scalar);

      double blend = 0.5 + 0.5 * scalar;
      if (blend < 0) blend = 0; if (blend > 1) blend = 1;

      unsigned char r = (1 - blend) * 60  + blend * 220;
      unsigned char g = (1 - blend) * 0   + blend * 20;
      unsigned char b = (1 - blend) * 100 + blend * 60;

      color->SetTuple3(i, r, g, b);
    }

    auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->SetPoints(points);

    std::vector<std::size_t> &tetraNodes = elementNodeTags[tetraIndex];
    int tetraCount = tetraNodes.size() / 4;

    for (int i = 0; i < tetraCount; i++) {
      auto tetra = vtkSmartPointer<vtkTetra>::New();
      for (int j = 0; j < 4; j++) {
        std::size_t gmshTag = tetraNodes[4 * i + j];
        tetra->GetPointIds()->SetId(j, tagToIndex[gmshTag]);
      }
      grid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
    }

    grid->GetPointData()->AddArray(velocity);
    grid->GetPointData()->SetActiveVectors("Velocity");
    grid->GetPointData()->AddArray(wave);
    grid->GetPointData()->AddArray(color);
    grid->GetPointData()->SetActiveScalars("Wave");

    auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(("output/frame_" + std::to_string(f) + ".vtu").c_str());
    writer->SetInputData(grid);
    writer->SetDataModeToAscii();
    writer->SetCompressor(nullptr);
    writer->Write();
  }

  gmsh::finalize();
  return 0;
}
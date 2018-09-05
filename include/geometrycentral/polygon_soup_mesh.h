#pragma once

#include <cstdlib>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <vector>

#include "geometrycentral/geometry.h"
#include <geometrycentral/utilities.h>
#include "geometrycentral/vector3.h"

namespace geometrycentral {

class PolygonSoupMesh {
 public:
  PolygonSoupMesh();
  PolygonSoupMesh(std::string meshFilename);
  PolygonSoupMesh(const std::vector<std::vector<size_t>>& polygons_,
                  const std::vector<Vector3>& vertexCoordinates_);

  // Mutate this mesh and by naively triangulating polygons
  void triangulate();

  // Mesh data
  std::vector<std::vector<size_t>> polygons;
  std::vector<std::vector<size_t>> uvPolygons;
  std::vector<Vector3> vertexCoordinates;
  std::vector<Vector2> uvCoordinates;

 private:
  void readMeshFromFile(std::string filename);
};

}  // namespace geometrycentral
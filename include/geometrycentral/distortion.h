#pragma once

#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/geometry.h"

using namespace geometrycentral;

class Distortion {
    public: 
        // computes (min, max, average) area distortion
        static Vector3 computeAreaScaling(HalfedgeMesh* mesh, Geometry<Euclidean>* geom);

        // computes (min, max, average) quasi conformal error
        static Vector3 computeQuasiConformalError(HalfedgeMesh* mesh, Geometry<Euclidean>* geom);

    private:
        static std::vector<double> distortion;
};
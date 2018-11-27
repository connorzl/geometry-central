#pragma once

#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/geometry.h"

using namespace geometrycentral;

class Distortion {
    public: 

        Distortion(HalfedgeMesh* mesh, Geometry<Euclidean>* geom);

        // computes (min, max, average) area distortion
        Vector3 computeAreaScaling();

        // computes (min, max, average) quasi conformal error
        Vector3 computeQuasiConformalError();

        bool computeGlobalOverlap();

        size_t computeTriangleFlips();

        float computeSeamLength();

        FaceData<double> areaDistortion;
        FaceData<double> angleDistortion;
        FaceData<double> trianglesFlipped;

    private:
        std::vector<double> distortion;
        HalfedgeMesh* mesh;
        Geometry<Euclidean>* geom;
};
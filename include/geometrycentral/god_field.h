#pragma once

#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/geometry.h"
#include "geometrycentral/operators.h"
#include "geometrycentral/linear_solvers.h"

using namespace geometrycentral;

class GodField {
    public:
        GodField() {};
        GodField(HalfedgeMesh *m, Geometry<Euclidean> *geom);

        void computeCrossField(EdgeData<double> &lengths, HalfedgeData<double> &angles);
        size_t computeSingularities(VertexData<int> &singularities);

        FaceData<std::complex<double>> field;

    private:
        HalfedgeMesh* mesh;
        Geometry<Euclidean> *geom;

        int nPowerIterations = 20;

        void setup();
        
        Eigen::SparseMatrix<std::complex<double>> assembleM();
        Eigen::SparseMatrix<std::complex<double>> assembleA();

        HalfedgeData<std::complex<double>> r;
        HalfedgeData<std::complex<double>> edgeVector;

        EdgeData<double> edgeLengths;
        HalfedgeData<double> heAngles;
        FaceData<double> faceAreas;
};
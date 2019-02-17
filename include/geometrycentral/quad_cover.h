#pragma once

#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/geometry.h"

using namespace geometrycentral;

class QuadCover {
    public:
        QuadCover(HalfedgeMesh* m, Geometry<Euclidean>* g);
        VertexData<std::complex<double>> computeCrossField();
        FaceData<int> computeSingularities();

    private:
        HalfedgeMesh* mesh;
        Geometry<Euclidean>* geom;

        int nPowerIterations = 20;
        double eps = std::pow(10.,-8.);

        // quantities
        HalfedgeData<std::complex<double>> r;
        VertexData<std::complex<double>> field;
        FaceData<int> singularities;

        // helpers
        HalfedgeData<std::complex<double>> setup();
        Eigen::SparseMatrix<std::complex<double>> assembleM();
        Eigen::SparseMatrix<std::complex<double>> assembleA();
        void computeSmoothestField(Eigen::SparseMatrix<std::complex<double>> M, Eigen::SparseMatrix<std::complex<double>> A);
};
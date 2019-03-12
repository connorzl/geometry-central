#pragma once

#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/geometry.h"
#include "geometrycentral/polygon_soup_mesh.h"

using namespace geometrycentral;

class QuadCover {
    public:
        QuadCover(HalfedgeMesh* m, Geometry<Euclidean>* g);
        VertexData<std::complex<double>> computeCrossField();
        FaceData<int> computeSingularities();
        void computeBranchCover();
        VertexData<double> computeOffset();

    private:
        HalfedgeMesh* mesh;
        Geometry<Euclidean>* geom;

        int nPowerIterations = 20;
        double eps = std::pow(10.,-8.);
        double n = 4;

        // quantities
        HalfedgeData<double> phi;
        HalfedgeData<std::complex<double>> r;
        VertexData<std::complex<double>> field;
        FaceData<int> singularities;
        HalfedgeData<int> branchCover;

        // helpers
        void setup();
        Eigen::SparseMatrix<std::complex<double>> assembleM();
        Eigen::SparseMatrix<std::complex<double>> assembleA();
        void computeSmoothestField(Eigen::SparseMatrix<std::complex<double>> M, Eigen::SparseMatrix<std::complex<double>> A);

        // visualization
        Eigen::SparseMatrix<std::complex<double>> buildLaplacian();
        void emitTriangles(VertexData<double> offset);
};
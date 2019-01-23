#pragma once

#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/geometry.h"
#include "geometrycentral/spectral_conformal.h"
#include <fstream> 

using namespace geometrycentral;

class ARAP {
    public:
        ARAP(HalfedgeMesh* m, Geometry<Euclidean>* g);
        void computeARAP();

    private:
        HalfedgeMesh* mesh;
        Geometry<Euclidean>* geom;
        VertexData<size_t> vertexIndices;
        HalfedgeData<Vector2> isoTriangleParam;
        VertexData<Vector2> uvCoords;

        // helpers
        void computeIsoTriangleParam();
        Eigen::SparseMatrix<std::complex<double>> createLaplaceMatrix();
        FaceData<Eigen::Matrix2d> computeLRotations(VertexData<Vector2> const &u);
        Eigen::MatrixXcd computebVector(FaceData<Eigen::Matrix2d> const &L);
        void writeToFile(std::ofstream &outfile);
        void normalize();
};
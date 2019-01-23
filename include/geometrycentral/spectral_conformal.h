#pragma once

#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/geometry.h"
#include <fstream> 

using namespace geometrycentral;

class SpectralConformal {
    public:
        SpectralConformal(HalfedgeMesh* m, Geometry<Euclidean>* g);
        VertexData<Vector2> computeSpectralConformal();

    private:
        HalfedgeMesh* mesh;
        Geometry<Euclidean>* geom;
        VertexData<size_t> vertexIndices;
        VertexData<Vector2> uvCoords;

        // helpers
        Eigen::SparseMatrix<std::complex<double>> createEDMatrix();
        Eigen::SparseMatrix<std::complex<double>> createAMatrix();
        double computeResidual(Eigen::SparseMatrix<std::complex<double>> EC, Eigen::MatrixXcd x);
        void writeToFile(std::ofstream &outfile);
        void normalize();
};
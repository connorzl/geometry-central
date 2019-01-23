#pragma once

#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/geometry.h"
#include <fstream> 

using namespace geometrycentral;

class LSCM {
    public:
        LSCM(HalfedgeMesh* m, Geometry<Euclidean>* g);
        VertexData<Vector2> computeLSCM();

    private:
        HalfedgeMesh* mesh;
        Geometry<Euclidean>* geom;
        VertexData<size_t> vertexIndices;
        VertexData<Vector2> uvCoords;

        // helpers
        void separateVertices(VertexPtr v1, VertexPtr v2);
        Eigen::SparseMatrix<std::complex<double>> createEDMatrix();
        Eigen::SparseMatrix<std::complex<double>> createAMatrix();
        void writeToFile(std::ofstream &outfile);
        void normalize();
};
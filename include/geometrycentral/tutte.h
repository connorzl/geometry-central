#pragma once

#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/geometry.h"
#include <assert.h> 
#include <fstream> 

using namespace geometrycentral;

class Tutte {
    public:
        Tutte(HalfedgeMesh* m, Geometry<Euclidean>* g);
        void computeTutteEmbedding();

    private:
        HalfedgeMesh* mesh;
        Geometry<Euclidean>* geom;
        VertexData<size_t> vertexIndices;
        VertexData<Vector2> uvCoords;

        // helpers
        void separateVertices();
        Eigen::MatrixXd mapToUnitCircle();
        Eigen::SparseMatrix<double> createKMatrix();
        void writeToFile(std::ofstream &outfile);
        void normalize();
};
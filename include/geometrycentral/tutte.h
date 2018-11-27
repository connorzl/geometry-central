#pragma once

#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/geometry.h"
#include <assert.h> 

using namespace geometrycentral;

class Tutte {
    public:
        Tutte(HalfedgeMesh* m, Geometry<Euclidean>* g);
        void computeTutteEmbedding();

    private:
        HalfedgeMesh* mesh;
        Geometry<Euclidean>* geom;

        // helpers
        void separateVertices(std::vector<VertexPtr> &interiorVertices, std::vector<VertexPtr> &boundaryVertices);
};
#include "geometrycentral/tutte.h"

using namespace geometrycentral;

Tutte::Tutte(HalfedgeMesh* m, Geometry<Euclidean>* g) : mesh(m), geom(g) {}

void Tutte::separateVertices(std::vector<VertexPtr> &interiorVertices, std::vector<VertexPtr> &boundaryVertices) {
    size_t numBoundary = 0;
    size_t numInterior = 0;

     // separate the vertices of the longest boundary loop from other boundary vertices
    size_t longestIndex = mesh->longestBoundaryLoop();
    for (size_t i = 0; i < mesh->nBoundaryLoops(); i++) {
        BoundaryPtr b = mesh->boundaryLoop(i);
        HalfedgePtr he = b.halfedge();
        HalfedgePtr curr = he;
        do {
            assert(!curr.isReal() && curr.vertex().isBoundary());
            
            VertexPtr v = curr.vertex();
            if (i == longestIndex) {
                boundaryVertices.push_back(v);
                numBoundary++;
            } else {
                interiorVertices.push_back(v);
                numInterior++;
            }
            curr = curr.next();
        } while(curr != he);
    }

    // now add in all the interior vertices
    size_t currVertex = 0;
    for (VertexPtr v : mesh->vertices()) {
        if (!v.isBoundary()) {
            interiorVertices.push_back(v);
            numInterior++;
        }
        currVertex++;
    }

    //std::cout<<boundaryVertices.size()<<std::endl;
    //std::cout<<interiorVertices.size()<<std::endl;
}
void Tutte::computeTutteEmbedding() {
    // Extract interior and exterior vertices
    std::vector<VertexPtr> interiorVertices;
    std::vector<VertexPtr> boundaryVertices;
    separateVertices(interiorVertices, boundaryVertices);

    // Map exterior vertices to unit circle

    // Build K matrix

    // Solve for interior vertex uv-coordinates
    
}
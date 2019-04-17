#include "geometrycentral/BHalfedgemesh.h"

std::vector<BFace> BranchCoverTopology::allFaces() {
    std::vector<BFace> BFaces; 
    for (FacePtr f : mesh->faces()) {
        for (int sheet = 0; sheet < 4; sheet++) {
            BFaces.push_back(BFace{ f, sheet, this });
        }
    }
    return BFaces;
}

std::vector<BVertex> BranchCoverTopology::allVertices() {
    std::vector<BVertex> BVertices;
    for (VertexPtr v : mesh->vertices()) {
        if (singular[v] != 0) {
            BVertices.push_back(BVertex{ v, singularSheet, this});
        } else {
            for (int sheet = 0; sheet < 4; sheet++) {
                BVertices.push_back(BVertex{ v, sheet, this });
            }
        }
    }
    return BVertices;
}

std::vector<BHalfedge> BranchCoverTopology::allHalfedges() {
    std::vector<BHalfedge> BHalfedges;
    for (HalfedgePtr he : mesh->allHalfedges()) {
        for (int sheet = 0; sheet < 4; sheet++) {
            BHalfedges.push_back(BHalfedge{ he, sheet, this });
        }
    }
    return BHalfedges;
}

std::vector<BEdge> BranchCoverTopology::allEdges() {
    std::vector<BEdge> BEdges;
    for (EdgePtr e : mesh->edges()) {
        for (int sheet = 0; sheet < 4; sheet++) {
            BEdges.push_back(BEdge{ e, sheet, this });
        }
    }
    return BEdges;
}

bool BranchCoverTopology::validateConnectivity() {
    std::cout << "Validating Connectivity...";

    //TODO 

    std::cout << "Done!" << std::endl;
    return true;
}
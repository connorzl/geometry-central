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
    std::cout << "Validating BHalfedgeMesh Connectivity...";

    // check edge and twin sanity
    for (BHalfedge he : allHalfedges()) {
        //std::cout << he.sheet << "," << he.twin().twin().sheet << std::endl;
        if (he != he.twin().twin()) throw std::logic_error("twins not reflective");
        if (he == he.twin()) throw std::logic_error("self-twin"); 
        if (he != he.edge().halfedge() && he.twin() != he.edge().halfedge()) throw std::logic_error("edge.halfedge doesn't match halfedge.edge");
        if (he.edge() != he.twin().edge()) throw std::logic_error ("he.edge and he.twin.edge are not the same");
    }

    for (BEdge e : allEdges()) {
        if (e.halfedge().edge() != e) throw std::logic_error("e.halfedge.edge != e");
    }

    // check face and next sanity
    size_t numHalfedges = allHalfedges().size();
    for (BFace f : allFaces()) {
        BHalfedge currHe = f.halfedge();
        BHalfedge firstHe = f.halfedge();
        size_t count = 0;
        do {
            if (currHe.face() != f) throw std::logic_error("face.halfedge doesn't match halfedge.face");
            currHe = currHe.next();
            count++;
            if (count > numHalfedges) throw std::logic_error ("next forms non-face loop");
        } while (currHe != firstHe);

        if (count < 3) throw std::logic_error("face of degree < 2");
    }

    for (BHalfedge he : allHalfedges()) {
        BHalfedge currHe = he;
        BHalfedge firstHe = he;
        size_t count = 0;
        do {
            if (currHe.face() != he.face()) throw std::logic_error("he.next.**.face doesn't match he.face");
            currHe = currHe.next();
            count ++;
            if (count > numHalfedges) throw std::logic_error("next forms non-face loop");
        } while (currHe != firstHe);

        if (he.vertex() != he.twin().next().vertex()) throw std::logic_error("halfedge vertices don't match");
    }

    // check vertex orbit sanity
    for (BVertex v : allVertices()) {
        BHalfedge currHe = v.halfedge();
        BHalfedge firstHe = v.halfedge();
        size_t count = 0;
        do {
            if (currHe.vertex() != v) throw std::logic_error("vertex.halfedge doesn't match halfedge.vertex");
            currHe = currHe.twin().next();
            count++;
            if (count > numHalfedges) throw std::logic_error("twin->next forms non-vertex loop");
        } while (currHe != firstHe);
    }

    std::cout << "Done!" << std::endl;
    return true;
}
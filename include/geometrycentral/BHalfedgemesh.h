
#pragma once

#include "geometrycentral/halfedge_mesh.h"

using namespace geometrycentral;

class BranchCoverTopology {
public:
    BranchCoverTopology(HalfedgeMesh* _mesh) : mesh(_mesh) {}

    // data
    HalfedgeMesh* mesh;
    HalfedgeData<int> sheetInterchange;
    VertexData<bool> isSingular;

    std::vector<BFace> allFaces();
    std::vector<BVertex> allVertices();
    std::vector<BEdge> allEdges();
    // ... more

    bool validateConnectivity();
};

struct BHalfedge {
    HalfedgePtr he;
    int sheet;
    BranchCoverTopology* BC;

    BHalfedge next();
    BHalfedge twin();
    BVertex vertex();
    BEdge edge();
    BFace face();
};

struct BVertex {
    VertexPtr v;
    int sheet;
    BranchCoverTopology* BC;

    BHalfedge halfedge(); 
};

struct BFace {
    FacePtr f;
    int sheet;
    BranchCoverTopology* BC;

    BHalfedge halfedge();
};

struct BEdge {
    EdgePtr e;
    int sheet;
    BranchCoverTopology* BC;

    BHalfedge halfedge();    
};


// equality operators

// BFace functions
inline BHalfedge BFace::halfedge() {
    return BHalfedge{ f.halfedge(), sheet, BC };
}

// BEdge functions
inline BHalfedge BEdge::halfedge() {
    return BHalfedge{ e.halfedge(), sheet, BC };
}

// BVertex functions
inline BHalfedge BVertex::halfedge() {
    if (BC->isSingular[v]) return BHalfedge{ v.halfedge(), 0, BC };
    return BHalfedge{ v.halfedge(), sheet, BC };
}

// BHalfedge functions
inline BHalfedge BHalfedge::next() {
    return BHalfedge{ he.next(), sheet, BC };
}

inline BHalfedge BHalfedge::twin() {
    return BHalfedge{ he.twin(), (sheet + BC->sheetInterchange[he]) % 4, BC };
}

inline BVertex BHalfedge::vertex() {
    // make sure that each outgoing BHalfedge has the same BVertex
    if (he == he.vertex().halfedge()) {
        return BVertex{ he.vertex(), sheet, BC };
    }
    return twin().next().vertex();
}

inline BEdge BHalfedge::edge() {
    // make sure that he.edge() == he.twin().edge()
    if (he == he.edge().halfedge()) {
        return BEdge{ he.edge(), sheet, BC };
    } else {
        return twin().edge();
    }
}

inline BFace BHalfedge::face() {
    return BFace{ he.face(), sheet, BC };
}
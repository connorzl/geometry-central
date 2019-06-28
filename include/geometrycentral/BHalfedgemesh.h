
#pragma once

#include "geometrycentral/halfedge_mesh.h"
#include <map>
#include <queue>

using namespace geometrycentral;

struct BFace;
struct BEdge;
struct BHalfedge;
struct BVertex;

class BranchCoverTopology {
public:
    BranchCoverTopology(){}
    BranchCoverTopology(HalfedgeMesh* _mesh, HalfedgeData<int> _eta, VertexData<int> _isSingular) : mesh(_mesh), eta(_eta), singularities(_isSingular) {}
    BranchCoverTopology(HalfedgeMesh* _mesh, VertexData<int> _singularities);

    // data
    HalfedgeMesh* mesh;
    HalfedgeData<int> eta;
    VertexData<int> singularities;
    int singularSheet = 0;

    std::vector<BFace> allFaces();
    std::vector<BVertex> allVertices();
    std::vector<BEdge> allEdges();
    std::vector<BHalfedge> allHalfedges();

    void computeEta();
    void validateEta(HalfedgeData<std::complex<double>> edgeVector);
    bool validateConnectivity();

    // computing sheet interchange function eta directly from singular points
    VertexPtr buildPrimalSpanningTree();
    std::map<HalfedgePtr,bool> inPrimalTree;
    std::map<VertexPtr,VertexPtr> vertexParent;
};

struct BFace {
    FacePtr f;
    int sheet;
    BranchCoverTopology* BC;

    BHalfedge halfedge();

    bool operator==(const BFace& other) const;
    bool operator!=(const BFace& other) const;
    bool operator<(const BFace& other) const;
};

struct BEdge {
    EdgePtr e;
    int sheet;
    BranchCoverTopology* BC;

    BHalfedge halfedge();   
    
    bool operator==(const BEdge& other) const;
    bool operator!=(const BEdge& other) const; 
};

struct BVertex {
    VertexPtr v;
    int sheet;
    BranchCoverTopology* BC;

    BHalfedge halfedge(); 
    std::vector<BFace> adjacentBFaces();

    bool operator==(const BVertex& other) const;
    bool operator!=(const BVertex& other) const;
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

    bool operator==(const BHalfedge& other) const;
    bool operator!=(const BHalfedge& other) const;
};

// BFace functions
inline BHalfedge BFace::halfedge() {
    return BHalfedge{ f.halfedge(), sheet, BC };
}

inline bool BFace::operator==(const BFace& other) const {
    return (f == other.f && sheet == other.sheet);
}

inline bool BFace::operator!=(const BFace& other) const {
    return !(*this == other);
}

inline bool BFace::operator<(const BFace& other) const {
    if (f == other.f) {
        return sheet < other.sheet;
    }
    return (f < other.f);
}

// BEdge functions
inline BHalfedge BEdge::halfedge() {
    return BHalfedge{ e.halfedge(), sheet, BC };
}

inline bool BEdge::operator==(const BEdge& other) const {
    return (e == other.e && sheet == other.sheet);
}

inline bool BEdge::operator!=(const BEdge& other) const {
    return !(*this == other);
}

// BVertex functions
inline BHalfedge BVertex::halfedge() {
    if (BC->singularities[v] != 0) return BHalfedge{ v.halfedge(), BC->singularSheet, BC };
    return BHalfedge{ v.halfedge(), sheet, BC };
}

inline std::vector<BFace> BVertex::adjacentBFaces() {
    std::vector<BFace> BFaces;

    // TODO: might need to fix this for meshes with boundary
    BHalfedge BHe_orig = halfedge();
    BHalfedge BHe_curr = BHe_orig;
    do {
        BFaces.push_back(BHe_curr.face());
        BHe_curr = BHe_curr.next().next().twin();
    } while (BHe_curr != BHe_orig);

    return BFaces;
}

inline bool BVertex::operator==(const BVertex& other) const {
    if (BC->singularities[v] == 0 && BC->singularities[other.v] == 0) {
        return ((v == other.v) && (sheet == other.sheet));
    } else if (BC->singularities[v] != 0 && BC->singularities[other.v] != 0) {
        if (sheet != other.sheet) throw std::logic_error("sheets of singular vertices do not match");
        return (v == other.v);
    } else {
        return false;
    }
}

inline bool BVertex::operator!=(const BVertex& other) const {
    return !(*this == other);
}

// BHalfedge functions
inline BHalfedge BHalfedge::next() {
    return BHalfedge{ he.next(), sheet, BC };
}

inline BHalfedge BHalfedge::twin() {
    if (!he.isReal() && BC->eta[he] != 0) {
        throw std::logic_error("Boundary halfedge eta is non-zero");
    }
    return BHalfedge{ he.twin(), (sheet + BC->eta[he]) % 4, BC };
}

inline BVertex BHalfedge::vertex() {    
    /*
    if (!he.isReal() && BC->eta[he] != 0) {
        throw std::logic_error("Boundary halfedge eta is non-zero");
    }
    if (he == he.vertex().halfedge()) {
        int s = (BC->singular[he.vertex()] != 0) ? BC->singularSheet : sheet;
        return BVertex{ he.vertex(), s, BC };
    }
    return twin().next().vertex();
    */
    // make sure that each outgoing BHalfedge has the same BVertex
    if (!he.isReal()) throw std::logic_error("imaginary BHalfedge being used for vertex");

    BHalfedge BHe_orig = BHalfedge(*this);
    BHalfedge BHe = BHe_orig;
    while (true) {
        if (!BHe.he.isReal()) break;
        if (BHe.he == BHe.he.vertex().halfedge()) {
            int s = (BC->singularities[BHe.he.vertex()] != 0) ? BC->singularSheet : BHe.sheet;
            return BVertex{ BHe.he.vertex(), s, BC };
        }

        BHe = BHe.twin().next();
    }

    BHe = BHe_orig;
    while (true) {
        if (!BHe.he.isReal()) break;

        if (BHe.he == BHe.he.vertex().halfedge()) {
            int s = (BC->singularities[BHe.he.vertex()] != 0) ? BC->singularSheet : BHe.sheet;
            return BVertex{ BHe.he.vertex(), s, BC };
        }

        BHe = BHe.next().next().twin();
    }
    throw std::logic_error("BHalfedge.vertex() code is broken rip");
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
    if (!he.isReal()) throw std::logic_error("imaginary BHalfedge being used for face");
    return BFace{ he.face(), sheet, BC };
}

inline bool BHalfedge::operator==(const BHalfedge& other) const {
    return (he == other.he && sheet == other.sheet);
}

inline bool BHalfedge::operator!=(const BHalfedge& other) const {
    return !(*this == other);
}
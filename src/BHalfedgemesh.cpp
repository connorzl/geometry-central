#include "geometrycentral/BHalfedgemesh.h"

BranchCoverTopology::BranchCoverTopology(HalfedgeMesh* m, VertexData<int> s) : mesh(m), singularities(s) {
    eta = HalfedgeData<int>(mesh,0);
    computeEta();
}

void BranchCoverTopology::computeEta() {
    VertexPtr root = buildPrimalSpanningTree();
   
    // count the number of children that each vertex has in the spanning tree
    VertexData<size_t> numChildren(mesh,0);
    for (VertexPtr v : mesh->vertices()) {
        size_t children = 0;
        for (HalfedgePtr he : v.outgoingHalfedges()) {
            if (inPrimalTree[he]) {
                children++;
            }
        }
        numChildren[v] = children;
    }

    // Queue of vertices that have all but 1 edge determined, starts with leaf nodes
    std::queue<VertexPtr> Q; 
    for (VertexPtr v : mesh->vertices()) {
        if (numChildren[v] == 0) {
            Q.push(v);
        }
    }

    size_t iVert = 0;
    while(!Q.empty()) {
        VertexPtr currVertex = Q.front();
        Q.pop();
        iVert++;

        int sheetSum = 0;
        for (HalfedgePtr he : currVertex.incomingHalfedges()) {
            sheetSum = (sheetSum + eta[he]) % 4;
        }

        if (currVertex == root) {
            if (sheetSum != singularities[root]) {
                std::cout << sheetSum << "," << singularities[root] << std::endl;
                throw std::logic_error("bad sheetSum around root vertex");
            } else {
                continue;
            }
        }

        for (HalfedgePtr he : currVertex.incomingHalfedges()) {
            if (inPrimalTree[he]) {
                if (he.vertex() != vertexParent[currVertex]) {
                    throw std::logic_error("tree broken");
                }

                int s = singularities[currVertex];
                if (s == 1 || s == -1 || s == 0) {
                    int diff = singularities[currVertex] - sheetSum;
                    while (diff < 0) {
                        diff += 4;
                    }

                    int offsetHe = diff % 4;
                    int offsetHeTwin = (4 - offsetHe) % 4;
                    
                    eta[he] = offsetHe;
                    eta[he.twin()] = offsetHeTwin;
                } else {
                    throw std::logic_error("only implemented for +1/-1 singularities");
                }
            }
        }
        // update parent
        VertexPtr parent = vertexParent[currVertex];
        numChildren[parent] -= 1;
        if (numChildren[parent] == 0) {
            Q.push(parent);
        }
    }

    if (iVert != mesh->nVertices()) {
        std::cout << iVert << "/" << mesh->nVertices() << std::endl;
        throw std::logic_error("did not visit all vertices");
    }
}

void BranchCoverTopology::validateEta(HalfedgeData<std::complex<double>> edgeVector) {
    for (VertexPtr v : mesh->vertices()) {
        if (v.isBoundary()) continue;
        int total = 0;
        for (HalfedgePtr he : v.incomingHalfedges()) {
            if (he.edge().isBoundary()) {
                eta[he] = 0;
                continue;
            }
            total = (total + eta[he]) % 4;
        }
       
        if (singularities[v] == 1 && total != 1) {
            std::cout << "difference at +1 singularity: " << total << std::endl;
        } else if (singularities[v] == -1 && total != 3) {
            std::cout << "difference at -1 singularity: " << total << std::endl;
        } else if (singularities[v] == 0 && total != 0) {
            std::cout << "difference at non-singularity: " << total << std::endl; 
        }
    }

    std::vector<BVertex> allVertices = this->allVertices();
    for (BVertex Bv : allVertices) {
        if (Bv.v.isBoundary()) continue;
        std::complex<double> v0(((double) rand() / (RAND_MAX)), ((double) rand() / (RAND_MAX)));
        std::complex<double> vCurr = v0;

        // rotate vector around fan
        BHalfedge BHe = Bv.halfedge();
        BHalfedge firstHe = BHe;
        do {
            std::complex<double> theta_ij = edgeVector[BHe.he];
            std::complex<double> theta_ji = edgeVector[BHe.he.twin()];
            std::complex<double> r_ij = -theta_ji / theta_ij;
            r_ij = r_ij / std::abs(r_ij);
            vCurr = r_ij * vCurr;
            BHe = BHe.twin().next();
        } while (BHe.he != firstHe.he);

        // check that the resulting vector makes sense
        std::complex<double> i(0,1);
        std::complex<double> rotPos90 = std::exp(i * M_PI_2);
        std::complex<double> rotNeg90 = std::exp(i * -M_PI_2);
        if (singularities[Bv.v] == 0 && std::norm(vCurr - v0) > 1e-8) {
            throw std::logic_error("wrong rotation around non-singular vertex");
        } else if (singularities[Bv.v] == 1 && std::norm(rotPos90 * vCurr - v0) > 1e-8) {
            throw std::logic_error("wrong rotation around +1 singular vertex");
        } else if (singularities[Bv.v] == -1 && std::norm(rotNeg90 * vCurr - v0) > 1e-8) {
            throw std::logic_error("wrong rotation around -1 singular vertex");
        }
    }
}

VertexPtr BranchCoverTopology::buildPrimalSpanningTree() {
    std::map<VertexPtr,bool> visited;
    for (HalfedgePtr he : mesh->realHalfedges()) {
        inPrimalTree[he] = false;
    }
    for (VertexPtr v : mesh->vertices()) {
        visited[v] = false;
    }
    
    VertexPtr root = mesh->vertex(0);
    visited[root] = true;
    vertexParent[root] = nullptr;
    size_t treeVertices = 0;

    std::queue<VertexPtr> Q; 
    Q.push(root);
    while(!Q.empty()) {
        VertexPtr currVertex = Q.front();
        Q.pop();

        // get neighbors			
        for (HalfedgePtr he : currVertex.outgoingHalfedges()) {
            VertexPtr v = he.next().vertex();
            if (!visited[v]) {
                Q.push(v);                       // add all unvisited neighbors into queue
                visited[v] = true;               // mark added neighbors as visited
                vertexParent[v] = currVertex;    // update the primal spanning tree
                inPrimalTree[he] = true;
                treeVertices++;
            }
        }
    }
    if( treeVertices != mesh->nVertices()-1 ) {
        throw std::logic_error("built primal spanning tree incorrectly");
    }
    return root;
}

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
        if (singularities[v] != 0) {
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
        if (he.he.edge().isBoundary()) continue;
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
        if (!he.he.isReal()) continue;
        BHalfedge currHe = he;
        BHalfedge firstHe = he;
        size_t count = 0;
        do {
            if (currHe.face() != he.face()) throw std::logic_error("he.next.**.face doesn't match he.face");
            currHe = currHe.next();
            count ++;
            if (count > numHalfedges) throw std::logic_error("next forms non-face loop");
        } while (currHe != firstHe);

        // the Bhalfedge.vertex() function does a twin().next() loop, which may hit boundary edges for vertices on boundary
        if (!he.twin().he.isReal()) continue;
        if (he.vertex() != he.twin().next().vertex()) {
            std::cout << "sheets: " << he.vertex().sheet << "," << he.twin().next().vertex().sheet << std::endl;
            throw std::logic_error("halfedge vertices don't match");
        } 
    }

    // check vertex orbit sanity
    for (BVertex v : allVertices()) {
        if (v.v.isBoundary()) continue;
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
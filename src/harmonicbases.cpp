#include "geometrycentral/harmonicbases.h"
#include <Eigen/SparseCholesky>
#include <Eigen/Dense>
#include <queue>

HarmonicBases::HarmonicBases(HalfedgeMesh* m, Geometry<Euclidean>* g) : mesh(m), geom(g), orientation(m) {
    for (EdgePtr e : mesh->edges()) {
        orientation[e.halfedge()] = 1;
        orientation[e.halfedge().twin()] = -1;
    }
}

void HarmonicBases::buildPrimalSpanningTree() {
    std::map<VertexPtr,bool> visited;
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
        for (VertexPtr v : currVertex.adjacentVertices()) {
            if (visited.find(v) == visited.end()) {
                Q.push(v);                       // add all unvisited neighbors into queue
                visited[v] = true;               // mark added neighbors as visited
                vertexParent[v] = currVertex;    // update the primal spanning tree
                treeVertices++;
            }
        }
    }
    assert( treeVertices == mesh->nVertices()-1 );
}

bool HarmonicBases::inPrimalSpanningTree(HalfedgePtr he) {
    VertexPtr v1 = he.vertex();
    VertexPtr v2 = he.twin().vertex();

    // check if (v2, v1) or (v1,v2) is in the primal spanning tree
    return ((vertexParent.find(v1) != vertexParent.end() && v2 == vertexParent[v1]) ||
            (vertexParent.find(v2) != vertexParent.end() && v1 == vertexParent[v2]));
}

void HarmonicBases::buildDualSpanningCotree() {
    std::map<FacePtr,bool> visited;
    FacePtr root = mesh->face(0);
	visited[root] = true;
    faceParent[root] = nullptr;

	std::queue<FacePtr> Q; 
    Q.push(root);
    size_t treeFaces = 0;
    while(!Q.empty()) {
        FacePtr currFace = Q.front();
        Q.pop();

        HalfedgePtr he = currFace.halfedge();
        do {
            FacePtr neighbor = he.twin().face();
            if (visited.find(neighbor) == visited.end() && !inPrimalSpanningTree(he)) {
                Q.push(neighbor);                       // add all unvisited neighbors into queue
                visited[neighbor] = true;               // mark added neighbors as visited
                faceParent[neighbor] = currFace;        // update the dual spanning co-tree
                treeFaces++;
            }
            he = he.next();
        } while (he != currFace.halfedge());
    }
    assert( treeFaces == mesh->nFaces()-1 );
}

bool HarmonicBases::inDualSpanningCotree(HalfedgePtr he) {
    FacePtr f1 = he.face();
    FacePtr f2 = he.twin().face();

    // check if (f2, f1) or (f1,f2) is in the primal spanning tree
    return ((faceParent.find(f1) != faceParent.end() && f2 == faceParent[f1]) ||
            (faceParent.find(f2) != faceParent.end() && f1 == faceParent[f2]));   
}

void HarmonicBases::treeCotree() {
    buildPrimalSpanningTree();
    buildDualSpanningCotree();

    // find generators
    for (EdgePtr e : mesh->edges()) {
        HalfedgePtr he = e.halfedge();

        // there should be exactly 2g of these halfedges
        if ((!inDualSpanningCotree(he) && (!inPrimalSpanningTree(he)))) {
            FacePtr f1 = he.face();
            FacePtr f2 = he.twin().face();
            std::vector<HalfedgePtr> loop_forward, loop_reverse;
            loop_forward.push_back(he);

            std::map<FacePtr,bool> visited;
            // follow one endpoint ("forward" along v1->v2->...->root)
            FacePtr currFace = f2;
            do {
                visited[currFace] = true;
                currFace = faceParent[currFace];
            } while (currFace != nullptr);
            
            // follow other endpoint ("reverse" along root->...->v1)
            currFace = f1;
            FacePtr sharedFace;
            while (faceParent[currFace] != nullptr) {
                FacePtr prevFace = faceParent[currFace];
                for (HalfedgePtr h : currFace.adjacentHalfedges()) {
                    if (h.twin().face() == prevFace)  {
                        loop_reverse.push_back(h.twin());   
                        break;
                    }
                }
                if (visited.find(prevFace) != visited.end()) {
                    sharedFace = prevFace;
                    break;
                }
                currFace = prevFace;
            }

            currFace = f2;
            while (currFace != sharedFace) {
                FacePtr nextFace = faceParent[currFace];
                for (HalfedgePtr h : currFace.adjacentHalfedges()) {
                    if (h.twin().face() == nextFace) {
                        loop_forward.push_back(h);
                        break;
                    }
                }
                currFace = nextFace;
            }		
            loop_forward.insert(loop_forward.end(), loop_reverse.begin(), loop_reverse.end());
            generators.push_back(loop_forward);
        }
    }
    assert( 2 - mesh->eulerCharacteristic() >= 0 );
    assert( generators.size() == size_t(2 - mesh->eulerCharacteristic()) );
}

Eigen::MatrixXd HarmonicBases::buildClosedPrimalOneForm(std::vector<HalfedgePtr> generator) {
    Eigen::MatrixXd omega = Eigen::MatrixXd::Zero(mesh->nEdges(),1);
    EdgeData<size_t> edgeIndices = mesh->getEdgeIndices();
    for (size_t i = 0; i < generator.size(); i++) {
        HalfedgePtr he = generator[i];
        size_t index = edgeIndices[he.edge()];
        omega(index,0) = orientation[he];
    }
    return omega;
}

Eigen::MatrixXd HarmonicBases::computeExactComponent(Eigen::MatrixXd omega) {
    Eigen::MatrixXd b = d0T * hodge1 * omega;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    Eigen::MatrixXd alpha = solver.solve(b);

    Eigen::MatrixXd dAlpha = d0 * alpha;
    return dAlpha;
}

// implement the harmonic bases next
std::vector<Eigen::MatrixXd> HarmonicBases::compute() {
    // find generators
    treeCotree();

    // compute Laplacian, d0T, and hodge1
    A = Operators::laplaceMatrix(mesh,geom);
    d0 = Operators::buildExteriorDerivative0Form(mesh,geom);
    d0T = d0.transpose();
    hodge1 = Operators::buildHodgeStar1Form(mesh,geom);

    // compute bases
    EdgeData<size_t> edgeIndices = mesh->getEdgeIndices();
    for (size_t i = 0; i < generators.size(); i++) {
		Eigen::MatrixXd omega = buildClosedPrimalOneForm(generators[i]);

        // sanity check
        for (FacePtr f : mesh->faces()) {
            int sum = 0;				
            for (HalfedgePtr he : f.adjacentHalfedges()) {
                size_t index = edgeIndices[he.edge()];
                sum += orientation[he] * omega(index,0);
            }
            if (sum != 0) {
                std::cout << "Generator " << i << " sum: " << sum << std::endl;
            }
        }
        Eigen::MatrixXd dAlpha = computeExactComponent(omega);
        Eigen::MatrixXd gamma = omega - dAlpha;
        bases.push_back(gamma);
    }
    return bases;
}

FaceData<std::complex<double>> HarmonicBases::visualize() {
    Eigen::MatrixXd omega = bases[0];
    FaceData<std::complex<double>> X(mesh);
    EdgeData<size_t> edgeIndices = mesh->getEdgeIndices();

    for (FacePtr f : mesh->faces()) {
        Eigen::MatrixXd A(3,2);
        Eigen::MatrixXd b(3,1);

        // Gather elements
        HalfedgePtr he = f.halfedge();
        double l_ij = geom->length(he.edge());
        double l_ik = geom->length(he.prev().edge());
        double theta_a = geom->angle(he.next());

        Vector2 i = Vector2{0,0};
        Vector2 j = Vector2{l_ij,0};
        Vector2 k = Vector2{cos(theta_a) * l_ik, sin(theta_a) * l_ik}; 

        Vector2 ij = j - i;
        Vector2 jk = k - j;
        Vector2 ki = i - k;

        A(0,0) = ij.x;
        A(0,1) = ij.y;
        A(1,0) = jk.x;
        A(1,1) = jk.y;
        A(2,0) = ki.x;
        A(2,1) = ki.y;

        b(0,0) = orientation[he       ] * omega(edgeIndices[he.edge()       ],0);
        b(1,0) = orientation[he.next()] * omega(edgeIndices[he.next().edge()],0);
        b(2,0) = orientation[he.prev()] * omega(edgeIndices[he.prev().edge()],0);

        Eigen::MatrixXd X_ijk = A.colPivHouseholderQr().solve(b);
        std::complex<double> res(-X_ijk(0,0), -X_ijk(1,0));
        if (std::abs(res) == 0) {
            X[f] = 0;
        } else {
            X[f] = res / std::abs(res);
        }
    } 
    return X;
}
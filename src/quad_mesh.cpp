#include "geometrycentral/quad_mesh.h"
#include <Eigen/SparseCholesky>
#include <queue>

// ONLY WORKS FOR MESHES WITHOUT BOUNDARY
QuadMesh::QuadMesh(HalfedgeMesh* m, Geometry<Euclidean>* g) : mesh(m), geom(g), theta(m), r(m), field(m), 
                                                                singularities(m), branchCover(m),
                                                                edgeLengthsCM(m), thetaCM(m), rCM(m), fieldCM(m) {
    assert(mesh->nBoundaryLoops() == 0);
}

void QuadMesh::setup() {
    for (FacePtr f : mesh->faces()) {
        // Gather elements
        HalfedgePtr he = f.halfedge();
        double l_jl = geom->length(he.edge());
        double l_ij = geom->length(he.prev().edge());
        double theta_ijl = geom->angle(he.next()); // radians

        // Place first vertex at (0,0)
        Vector2 j = Vector2{0,0};
        Vector2 l = Vector2{l_jl,0};
        Vector2 i = Vector2{cos(theta_ijl) * l_ij, sin(theta_ijl) * l_ij}; 

        Vector2 jl = l - j;
        Vector2 li = i - l;
        Vector2 ij = j - i;
        jl.normalize();
        li.normalize();
        ij.normalize();

        theta[he] = std::complex<double>(jl.x,jl.y);
        theta[he.next()] = std::complex<double>(li.x,li.y);
        theta[he.prev()] = std::complex<double>(ij.x,ij.y);
    }
    
    // Compute d_ij * r_ji = -d_ji
    for (VertexPtr v : mesh->vertices()) {
        for (HalfedgePtr he : v.outgoingHalfedges()) {
            std::complex<double> theta_ij = theta[he];
            std::complex<double> theta_ji = theta[he.twin()];
            r[he] = std::pow((-theta_ij / theta_ji), n);
        }
    }   
}

Eigen::SparseMatrix<std::complex<double>> QuadMesh::assembleM() {
    size_t n = mesh->nFaces();
    Eigen::SparseMatrix<std::complex<double>> M(n,n);
    std::vector<Eigen::Triplet<std::complex<double>>> triplets;

    FaceData<size_t> faceIndices = mesh->getFaceIndices();
    for (FacePtr f : mesh->faces()) {
        size_t i = faceIndices[f];
        triplets.push_back(Eigen::Triplet<std::complex<double>>(i, i, geom->area(f)));
    }
    M.setFromTriplets(triplets.begin(),triplets.end());
    return M;
}

Eigen::SparseMatrix<std::complex<double>> QuadMesh::assembleA() {
    size_t n = mesh->nFaces();
    Eigen::SparseMatrix<std::complex<double>> A(n,n);
    std::vector<Eigen::Triplet<std::complex<double>>> triplets;

    FaceData<size_t> faceIndices = mesh->getFaceIndices();
    for (FacePtr f : mesh->faces()) {
        size_t i = faceIndices[f];

        std::complex<double> r_ij = r[f.halfedge()];
        std::complex<double> r_jk = r[f.halfedge().next()];
        std::complex<double> r_ki = r[f.halfedge().prev()];

        size_t i_ij = faceIndices[f.halfedge().twin().face()];
        size_t i_jk = faceIndices[f.halfedge().next().twin().face()];
        size_t i_ki = faceIndices[f.halfedge().prev().twin().face()];

        triplets.push_back(Eigen::Triplet<std::complex<double>>(i, i, 3));
        triplets.push_back(Eigen::Triplet<std::complex<double>>(i,i_ij,-r_ij));
        triplets.push_back(Eigen::Triplet<std::complex<double>>(i,i_jk,-r_jk));
        triplets.push_back(Eigen::Triplet<std::complex<double>>(i,i_ki,-r_ki));
    }
    
    A.setFromTriplets(triplets.begin(),triplets.end());
    return A;
}

void QuadMesh::computeSmoothestField(Eigen::SparseMatrix<std::complex<double>> M, Eigen::SparseMatrix<std::complex<double>> A) {
    // LL^T <- Cholesky(A)
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>> solver;
    solver.compute(A);

    // u <- UniformRand(-1,1)
    Eigen::MatrixXcd u = Eigen::MatrixXcd::Random(mesh->nFaces(),1);
    Eigen::MatrixXcd x;

    // inverse power iteration to find eigenvector belonging to the smallest eigenvalue
    for (int i = 0; i < nPowerIterations; i++) {
        x = solver.solve(M * u);
        std::complex<double> norm2 = (x.transpose() * M * x)(0,0);
        u = x / sqrt(norm2);
    }

    // map resulting vector to VertexData
    FaceData<size_t> faceIndices = mesh->getFaceIndices();
    for (FacePtr f : mesh->faces()) {
        std::complex<double> c = u(faceIndices[f],0);
        if (std::abs(c) == 0) {
            field[f] = 0;
        } else {
            field[f] = c / std::abs(c);   
        }
    }
} 

FaceData<std::complex<double>> QuadMesh::computeCrossField() {
    std::cout << "Computing Cross Field...";
    // Algorithm 1 : Setup
    setup();

    // Algorithm 2 : Smoothest Field
    Eigen::SparseMatrix<std::complex<double>> M = assembleM();
    Eigen::SparseMatrix<std::complex<double>> A = assembleA();
    A = A + eps * M;
    computeSmoothestField(M,A);

    std::cout << "Done!" << std::endl;
    return field;
}

VertexData<int> QuadMesh::computeSingularities() {
    std::cout << "Computing Singularities...";

    // finally, compute index for each vertex v
    int total = 0;
    for (VertexPtr v : mesh->vertices()) {
        double angleSum = 0;

        for (HalfedgePtr he : v.outgoingHalfedges()) {
            std::complex<double> u_i = field[he.face()];
            std::complex<double> u_j = field[he.twin().face()];
            std::complex<double> r_ji = r[he];
            angleSum += std::arg(u_i / (r_ji * u_j));
        }

        double phi = (angleSum + n*geom->angleDefect(v)) / (2.0 * M_PI);
        singularities[v] = std::round(phi);
        total += singularities[v];
    }
    std::cout << "Done! Singularities Index Sum: " << total << std::endl;
    return singularities;
}

HalfedgeData<int> QuadMesh::computeBranchCover() {
    std::cout<< "Computing Branch Cover...";
    std::complex<double> i(0, 1);
    for (VertexPtr v : mesh->vertices()) {
        int total = 0;
        for (HalfedgePtr he : v.outgoingHalfedges()) {
            std::complex<double> u_i = std::pow(field[he.face()], 1.0 / n);
            std::complex<double> u_j = std::pow(field[he.twin().face()], 1.0 / n);
            
            std::complex<double> theta_ij = theta[he];
            std::complex<double> theta_ji = theta[he.twin()];
            std::complex<double> r_ji = -theta_ij / theta_ji;
            double ang = std::arg(u_i / (r_ji * u_j));

            if (ang >= -M_PI_4 && ang < M_PI_4) {
                branchCover[he] = 0;
            } else if (ang >= M_PI_4 && ang < 3.0 * M_PI_4) {
                branchCover[he] = 1;
                total = (total + 1) % n;
            } else if ((ang >= 3.0 * M_PI_4 && ang <= PI) || 
                       (ang < -3.0 * M_PI_4 && ang >= -PI)) {
                branchCover[he] = 2;
                total = (total + 2) % n;
            } else {
                assert(ang >= -3.0 * M_PI_4 && ang < -M_PI_4);
                branchCover[he] = 3;
                total = (total + 3) % n;
            }
        }
        if (singularities[v] != 0 && total == 0) {
            std::cout << "difference at singularity: " << total << std::endl;
        } else if (singularities[v] == 0 && total != 0) {
            std::cout << "difference at non-singularity: " << total << std::endl; 
        }
    }
    std::cout << "Done!" << std::endl;
    return branchCover;
}

VertexData<double> QuadMesh::uniformize() {
    size_t n = mesh->nVertices();
    VertexData<size_t> vertexIndices = mesh->getVertexIndices();
    
    // Ax = b, where A = |V| x |V| laplacian, x = |V| x 1 vertex scaling, b = |V| x 1 curvature diff K - K*
    Eigen::MatrixXd KTarg(n,1);
    for (VertexPtr v : mesh->vertices()) {
        size_t index = vertexIndices[v];
        if (singularities[v] == 0) {
            KTarg(index,0) = 0;
        } else if (singularities[v] == 1) {
            KTarg(index,0) = M_PI_2;
        } else {
            assert(singularities[v] == -1);
            KTarg(index,0) = -M_PI_2;
        }
    } 

    EdgeData<double> l0(mesh);
    geom->getEdgeLengths(l0);
    geom->getEdgeLengths(edgeLengthsCM);

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    Eigen::MatrixXd u = Eigen::MatrixXd::Zero(n,1);
    Eigen::SparseMatrix<double> A;
    Eigen::MatrixXd K,x;
    do {
        A = Operators::intrinsicLaplaceMatrix(mesh,edgeLengthsCM);
        K = Operators::intrinsicCurvature(mesh,edgeLengthsCM);

        solver.compute(A);
        x = solver.solve(KTarg - K); // x.norm() becomes nan
        std::cout << x.norm() << std::endl;
        u = u + x;

        // update edge lengths
        for (EdgePtr e : mesh->edges()) {
            VertexPtr vi = e.halfedge().vertex();
            VertexPtr vj = e.halfedge().twin().vertex();

            double ui = u(vertexIndices[vi],0);
            double uj = u(vertexIndices[vj],0);

            double s = std::exp( (ui + uj) / 2 );
            edgeLengthsCM[e] = l0[e] * s;
        }

        // check to see if triangle inequality still holds
        bool triangleInequality = true;
        for (FacePtr f : mesh->faces()) {
            HalfedgePtr he = f.halfedge();
            double a = edgeLengthsCM[he.edge()];
            double b = edgeLengthsCM[he.next().edge()];
            double c = edgeLengthsCM[he.prev().edge()];
            if (a > b + c || b > a + c || c > a + b) {
                triangleInequality = false;
                std::cout << "TRIANGLE INEQUALITY VIOLATED!" << std::endl;
                break;
            }
        }    
        if (!triangleInequality) break;
    } while (x.norm() > 1e-5);

    // store curvatures for visualization
    VertexData<double> curvatures(mesh);
    for (VertexPtr v : mesh->vertices()) {
        curvatures[v] = K(vertexIndices[v],0);        
    }
    return curvatures;
}

void QuadMesh::setupCM() {
    HalfedgeData<double> angles = Operators::computeAngles(mesh, edgeLengthsCM); 
    
    for (FacePtr f : mesh->faces()) {
        // Gather elements
        HalfedgePtr he = f.halfedge();
        double l_jl = edgeLengthsCM[he.edge()];
        double l_ij = edgeLengthsCM[he.prev().edge()];
        double theta_ijl = angles[he.next()];

        // Place first vertex at (0,0)
        Vector2 j = Vector2{0,0};
        Vector2 l = Vector2{l_jl,0};
        Vector2 i = Vector2{cos(theta_ijl) * l_ij, sin(theta_ijl) * l_ij}; 

        Vector2 jl = l - j;
        Vector2 li = i - l;
        Vector2 ij = j - i;
        jl.normalize();
        li.normalize();
        ij.normalize();

        thetaCM[he] = std::complex<double>(jl.x,jl.y);
        thetaCM[he.next()] = std::complex<double>(li.x,li.y);
        thetaCM[he.prev()] = std::complex<double>(ij.x,ij.y);
    }
    
    // Compute d_ij * r_ji = -d_ji
    for (VertexPtr v : mesh->vertices()) {
        for (HalfedgePtr he : v.outgoingHalfedges()) {
            std::complex<double> theta_ij = theta[he];
            std::complex<double> theta_ji = theta[he.twin()];
            rCM[he] = -theta_ij / theta_ji;
        }
    }   
}

FaceData<std::complex<double>> QuadMesh::computeCrossFieldCM() {
    setupCM();

    // tile faces using BFS and the r_ij
    // keep in mind that r_ij should be used to transfrom a vector in the face associated with he_ji
    // so, to transport a vector x_ij across an edge, do x_ij * r[he_ij.twin()]
    std::map<FacePtr,bool> visited;
    FacePtr root = mesh->face(0);
    fieldCM[root] = std::complex<double>(0.5,0.5);

    std::queue<FacePtr> Q; 
    Q.push(root);
    while(!Q.empty()) {
        FacePtr currFace = Q.front();
        Q.pop();

        // get neighbors			
        for (HalfedgePtr he : currFace.adjacentHalfedges()) {
            FacePtr neighbor = he.twin().face();
            if (visited.find(neighbor) == visited.end()) {
                Q.push(neighbor);              // add all unvisited neighbors into queue
                visited[neighbor] = true;      // mark added neighbors as visited
                fieldCM[neighbor] = rCM[he.twin()] * fieldCM[currFace];
            }
        }
    }
    for (FacePtr f : mesh->faces()) {
        fieldCM[f] = std::pow(fieldCM[f],n);
    }
    return fieldCM;
}
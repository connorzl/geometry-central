#include "geometrycentral/quad_mesh.h"
#include <Eigen/SparseCholesky>

// ONLY WORKS FOR MESHES WITHOUT BOUNDARY
QuadMesh::QuadMesh(HalfedgeMesh* m, Geometry<Euclidean>* g) : mesh(m), geom(g), theta(m), r(m), field(m), 
                                                                singularities(m), branchCover(m) {
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

void QuadMesh::computeBranchCover() {
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
}
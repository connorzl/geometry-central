#include "geometrycentral/god_field.h"

GodField::GodField(HalfedgeMesh *m, Geometry<Euclidean>* g) : mesh(m), geom(g) {
    edgeVector = HalfedgeData<std::complex<double>>(mesh);
    r = HalfedgeData<std::complex<double>>(mesh);
    field = FaceData<std::complex<double>>(mesh);
    faceAreas = FaceData<double>(mesh);
} 

void GodField::computeCrossField(EdgeData<double> &lengths, HalfedgeData<double> &angles) {
    std::cout << "Computing Smoothest Cross Field" << std::endl;
    
    edgeLengths = lengths;
    heAngles = angles;
    faceAreas = Operators::computeAreas(mesh, edgeLengths);
    setup();

    // Get matrices
    Eigen::SparseMatrix<std::complex<double>> M = assembleM();
    Eigen::SparseMatrix<std::complex<double>> A = assembleA();
    A = A + 1e-8 * M;

    // u <- UniformRand(-1,1)
    Eigen::MatrixXcd u = Eigen::MatrixXcd::Random(mesh->nFaces(),1);
    Eigen::MatrixXcd x;

    // inverse power iteration to find eigenvector belonging to the smallest eigenvalue
    PositiveDefiniteSolver<std::complex<double>> s(A);
    for (int i = 0; i < nPowerIterations; i++) {
        Eigen::MatrixXcd rhs = M * u;
        x = s.solve(rhs);

        std::complex<double> norm2 = (x.transpose() * M * x)(0,0);
        u = x / sqrt(norm2);
    }

    // map resulting vector to VertexData
    FaceData<size_t> faceIndices = mesh->getFaceIndices();
    for (FacePtr f : mesh->faces()) {
        std::complex<double> c = u(faceIndices[f],0);
        std::complex<double> val;
        if (std::abs(c) == 0) {
            val = 0;
        } else {
            val = c / std::abs(c);   
        }

        field[f] = val;
    }
    std::cout << "Done!" << std::endl;
}

void GodField::setup() {
    for (FacePtr f : mesh->faces()) {
        // Gather elements
        HalfedgePtr he = f.halfedge();
        double l_jl = edgeLengths[he.edge()];
        double l_ij = edgeLengths[he.prev().edge()];
        double theta_ijl = heAngles[he.next()]; // radians

        // Place first vertex at (0,0)
        Vector2 j = Vector2{0,0};
        Vector2 l = Vector2{l_jl,0};
        Vector2 i = Vector2{cos(theta_ijl) * l_ij, sin(theta_ijl) * l_ij}; 

        Vector2 jl = l - j;
        Vector2 li = i - l;
        Vector2 ij = j - i;

        edgeVector[he] = std::complex<double>(jl.x,jl.y);
        edgeVector[he.next()] = std::complex<double>(li.x,li.y);
        edgeVector[he.prev()] = std::complex<double>(ij.x,ij.y);
    }
    
    // Compute d_ij * r_ji = -d_ji
    for (VertexPtr v : mesh->vertices()) {
        for (HalfedgePtr he : v.outgoingHalfedges()) {
            if (he.edge().isBoundary()) {
                continue;
            }
            std::complex<double> theta_ij = edgeVector[he];
            std::complex<double> theta_ji = edgeVector[he.twin()];
            r[he] = std::pow((-theta_ij / theta_ji), 4);
            r[he] = r[he] / std::abs(r[he]);
        }
    }   
}

Eigen::SparseMatrix<std::complex<double>> GodField::assembleM() {
    size_t n = mesh->nFaces();
    Eigen::SparseMatrix<std::complex<double>> M(n,n);
    std::vector<Eigen::Triplet<std::complex<double>>> triplets;

    FaceData<size_t> faceIndices = mesh->getFaceIndices();
    for (FacePtr f : mesh->faces()) {
        size_t i = faceIndices[f];
        triplets.push_back(Eigen::Triplet<std::complex<double>>(i, i, faceAreas[f]));
    }
    M.setFromTriplets(triplets.begin(),triplets.end());
    return M;
}

Eigen::SparseMatrix<std::complex<double>> GodField::assembleA() {
    size_t n = mesh->nFaces();
    Eigen::SparseMatrix<std::complex<double>> A(n,n);
    std::vector<Eigen::Triplet<std::complex<double>>> triplets;

    FaceData<size_t> faceIndices = mesh->getFaceIndices();
    for (FacePtr f : mesh->faces()) {
        size_t i = faceIndices[f];
        HalfedgePtr he_ij = f.halfedge();

        std::complex<double> r_ij, r_jk, r_ki;
        r_ij = r[he_ij];
        r_jk = r[he_ij.next()];
        r_ki = r[he_ij.prev()];

        double sumWeight = 0;
        double weight;
        HalfedgePtr he = f.halfedge();
    
        if (!he_ij.edge().isBoundary()) {
            weight = 1; // edgeLengths[he_ij.edge()];
            sumWeight += weight;
            size_t i_ij = faceIndices[he_ij.twin().face()];
            triplets.push_back(Eigen::Triplet<std::complex<double>>(i,i_ij,-weight * r_ij));
        }
        if (!he_ij.next().edge().isBoundary()) {
            weight = 1;//edgeLengths[he_ij.next().edge()];
            sumWeight += weight;
            size_t i_jk = faceIndices[he_ij.next().twin().face()];
            triplets.push_back(Eigen::Triplet<std::complex<double>>(i,i_jk,-weight * r_jk));
        }
        if (!he_ij.prev().edge().isBoundary()) {
            weight = 1;//edgeLengths[he_ij.prev().edge()];
            sumWeight += weight;
            size_t i_ki = faceIndices[he_ij.prev().twin().face()];
            triplets.push_back(Eigen::Triplet<std::complex<double>>(i,i_ki,-weight * r_ki));   
        }
        triplets.push_back(Eigen::Triplet<std::complex<double>>(i, i, sumWeight));
    }
    A.setFromTriplets(triplets.begin(),triplets.end());
    return A;
}

size_t GodField::computeSingularities(VertexData<int> &singularities) {
    std::cout << "Computing Singularities...";

    // compute index for each vertex v
    int total = 0;
    size_t numSingularities = 0;
    for (VertexPtr v : mesh->vertices()) {
        double angleSum = 0;

        if (v.isBoundary()) {
            singularities[v] = 0;
            continue;
        } 

        for (HalfedgePtr he : v.outgoingHalfedges()) {
            std::complex<double> u_i = field[he.face()];
            std::complex<double> u_j = field[he.twin().face()];
            std::complex<double> r_ji = r[he];
            angleSum += std::arg(u_i / (r_ji * u_j));
        }

        double phi = (angleSum + 4*geom->angleDefect(v)) / (2.0 * M_PI);
        singularities[v] = std::round(phi);
        if (singularities[v] != 0) {
            numSingularities++;
        }
        total += singularities[v];
    }
    std::cout << "Done! Singularities Index Sum: " << total << std::endl;
    return numSingularities;
}
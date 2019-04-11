#include "geometrycentral/quad_cover.h"
#include <Eigen/SparseCholesky>

// ONLY WORKS FOR MESHES WITHOUT BOUNDARY
QuadCover::QuadCover(HalfedgeMesh* m, Geometry<Euclidean>* g) : mesh(m), geom(g), phi(m), r(m), field(m), 
                                                                singularities(m), branchCover(m) {
    assert(mesh->nBoundaryLoops() == 0);
}

void QuadCover::setup() {
    VertexData<double> s(mesh);
    // Compute s_i at each vertex
    for (VertexPtr v : mesh->vertices()) {
        if (v.isBoundary()) {
            s[v] = 1.0;
        } else {
            double sum = 0;
            for (HalfedgePtr he : v.outgoingHalfedges()) {
                sum += geom->angle(he.next());
            }
            s[v] = 2*M_PI / sum;
        }
    }
    
    // Compute transport at edges r_ij <- e^ip_ij
    for (VertexPtr v : mesh->vertices()) {
        HalfedgePtr he = v.halfedge();
        double angle = 0;
        double s_i = s[v];
        do {
            phi[he] = angle;
            angle += s_i * geom->angle(he.next());
            he = he.next().next().twin();
        } while (he != v.halfedge());
    }

    // Compute r_ij
    std::complex<double> i(0, 1);
    for (VertexPtr v : mesh->vertices()) {
        for (HalfedgePtr he : v.outgoingHalfedges()) {
            double theta_ij = phi[he];
            double theta_ji = phi[he.twin()] + M_PI;
            double rho_ij = theta_ij - theta_ji;
            r[he] = std::exp(i * n * rho_ij);
        }
    }   
}

Eigen::SparseMatrix<std::complex<double>> QuadCover::assembleM() {
    size_t n = mesh->nVertices();
    Eigen::SparseMatrix<std::complex<double>> M(n,n);
    std::vector<Eigen::Triplet<std::complex<double>>> triplets;

    VertexData<size_t> vertexIndices = mesh->getVertexIndices();
    for (FacePtr f : mesh->faces()) {
        HalfedgePtr he_ij = f.halfedge();
        size_t i = vertexIndices[he_ij.vertex()];
        size_t j = vertexIndices[he_ij.next().vertex()];
        size_t k = vertexIndices[he_ij.prev().vertex()];

        double area = geom->area(f);
        triplets.push_back(Eigen::Triplet<std::complex<double>>(i, i, area/3.));
        triplets.push_back(Eigen::Triplet<std::complex<double>>(j, j, area/3.));
        triplets.push_back(Eigen::Triplet<std::complex<double>>(k, k, area/3.));
    }
    M.setFromTriplets(triplets.begin(),triplets.end());
    return M;
}

Eigen::SparseMatrix<std::complex<double>> QuadCover::assembleA() {
    size_t n = mesh->nVertices();
    Eigen::SparseMatrix<std::complex<double>> A(n,n);
    std::vector<Eigen::Triplet<std::complex<double>>> triplets;

    VertexData<size_t> vertexIndices = mesh->getVertexIndices();
    for (FacePtr f : mesh->faces()) {
        HalfedgePtr he_ij = f.halfedge();
        HalfedgePtr he_jk = f.halfedge().next();
        HalfedgePtr he_ki = f.halfedge().prev();

        size_t i = vertexIndices[he_ij.vertex()];
        size_t j = vertexIndices[he_jk.vertex()];
        size_t k = vertexIndices[he_ki.vertex()];

        double a = geom->cotan(he_jk);
        double b = geom->cotan(he_ki);
        double c = geom->cotan(he_ij);

        std::complex<double> r_ij = r[he_ij];
        std::complex<double> r_ji = r[he_ij.twin()];
        std::complex<double> r_jk = r[he_jk];
        std::complex<double> r_kj = r[he_jk.twin()];
        std::complex<double> r_ki = r[he_ki];
        std::complex<double> r_ik = r[he_ki.twin()];

        // row i
        triplets.push_back(Eigen::Triplet<std::complex<double>>(i,i,b + c));
        triplets.push_back(Eigen::Triplet<std::complex<double>>(i,j,-c * r_ij));
        triplets.push_back(Eigen::Triplet<std::complex<double>>(i,k,-b * r_ik));

        // row j
        triplets.push_back(Eigen::Triplet<std::complex<double>>(j,i,-c * r_ji));
        triplets.push_back(Eigen::Triplet<std::complex<double>>(j,j,c + a));
        triplets.push_back(Eigen::Triplet<std::complex<double>>(j,k,-a * r_jk));

        // row k
        triplets.push_back(Eigen::Triplet<std::complex<double>>(k,i,-b * r_ki));
        triplets.push_back(Eigen::Triplet<std::complex<double>>(k,j,-a * r_kj));
        triplets.push_back(Eigen::Triplet<std::complex<double>>(k,k,a + b));
    }
    A.setFromTriplets(triplets.begin(),triplets.end());
    return A;
}

void QuadCover::computeSmoothestField(Eigen::SparseMatrix<std::complex<double>> M, Eigen::SparseMatrix<std::complex<double>> A) {
    // LL^T <- Cholesky(A)
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>> solver;
    solver.compute(A);

    // u <- UniformRand(-1,1)
    Eigen::MatrixXcd u = Eigen::MatrixXcd::Random(mesh->nVertices(),1);
    Eigen::MatrixXcd x;

    // inverse power iteration to find eigenvector belonging to the smallest eigenvalue
    for (int i = 0; i < nPowerIterations; i++) {
        x = solver.solve(M * u);
        std::complex<double> norm2 = (x.transpose() * M * x)(0,0);
        u = x / sqrt(norm2);
    }

    // map resulting vector to VertexData
    VertexData<size_t> vertexIndices = mesh->getVertexIndices();
    for (VertexPtr v : mesh->vertices()) {
        std::complex<double> c = u(vertexIndices[v],0);
        if (std::abs(c) == 0) {
            field[v] = 0;
        } else {
            field[v] = c / std::abs(c);   
        }
    }
} 

VertexData<std::complex<double>> QuadCover::computeCrossField() {
    std::cout << "Computing Cross Field!" << std::endl;
    // Algorithm 1 : Setup
    setup();

    // Algorithm 2 : Smoothest Field
    Eigen::SparseMatrix<std::complex<double>> M = assembleM();
    Eigen::SparseMatrix<std::complex<double>> A = assembleA();
    A = A + eps * M;
    computeSmoothestField(M,A);

    std::cout << "Done computing smoothest field!" << std::endl;
    return field;
}

FaceData<int> QuadCover::computeSingularities() {
    std::cout << "Computing Singularities!" << std::endl;
    // first, compute Omega_ijk <- arg(r_ij r_jk r_ki)
    FaceData<double> Omega(mesh);
    for (FacePtr f : mesh->faces()) {
        std::complex<double> r_ij = r[f.halfedge()];
        std::complex<double> r_jk = r[f.halfedge().next()];
        std::complex<double> r_ki = r[f.halfedge().prev()];
        Omega[f] = std::arg(r_ij * r_jk * r_ki);
    }

    // next, compute w_ij for each e_ij, such that u_j = e^iw_ij * r_ij * u_i 
    // w_ij = arg(u_j / (r_ij * u_i))
    HalfedgeData<double> w(mesh);
    for (HalfedgePtr he : mesh->allHalfedges()) {
        std::complex<double> u_i = field[he.vertex()];
        std::complex<double> u_j = field[he.twin().vertex()];        
        std::complex<double> r_ij = r[he];
        w[he] = std::arg(u_j * r_ij / u_i);
    }

    // finally, compute index for each triangle t
    // (1/2pi) * (w_ij + w_jk + w_ki + Omega_ijk)
    for (FacePtr f : mesh->faces()) {
        double w_ij = w[f.halfedge()];
        double w_jk = w[f.halfedge().next()];
        double w_ki = w[f.halfedge().prev()];
        double Omega_ijk = Omega[f];
        double phi = (w_ij + w_jk + w_ki - Omega_ijk) / (2.0 * M_PI);
        singularities[f] = std::round(phi);
    }
    return singularities;
}

void QuadCover::computeBranchCover() {
    std::cout<< "Computing Branch Cover!" << std::endl;
    std::complex<double> i(0, 1);
    for (FacePtr f : mesh->faces()) {
        int total = 0;
        for (HalfedgePtr he_ij : f.adjacentHalfedges()) {
            HalfedgePtr he_ji = he_ij.twin();

            std::complex<double> f_ij = std::pow(field[he_ij.vertex()], 1.0 / n);
            std::complex<double> f_ji = std::pow(field[he_ji.vertex()], 1.0 / n);
            
            // we need to recompute r_ij here without raising to the nth power, 
            // as raising to the nth power and then taking the nth root is not always an identity operation
            double theta_ij = phi[he_ij];
            double theta_ji = phi[he_ij.twin()] + M_PI;
            double rho_ij = theta_ji - theta_ij;
            std::complex<double> r_ij = std::exp(i * rho_ij);
            std::complex<double> s_ij = f_ji / (f_ij * r_ij); 
            double ang = std::arg(s_ij);
            
            if (ang >= -M_PI_4 && ang < M_PI_4) {
                branchCover[he_ij] = 0;
            } else if (ang >= M_PI_4 && ang < 3.0 * M_PI_4) {
                branchCover[he_ij] = 1;
                total = (total + 1) % 4;
            } else if ((ang >= 3.0 * M_PI_4 && ang <= PI) || 
                       (ang < -3.0 * M_PI_4 && ang >= -PI)) {
                branchCover[he_ij] = 2;
                total = (total + 2) % 4;
            } else {
                assert(ang >= -3.0 * M_PI_4 && ang < -M_PI_4);
                branchCover[he_ij] = 3;
                total = (total + 3) % 4;
            }
            /*          
            if ( ang >= -M_PI_2 && ang < M_PI_2 ) {
                branchCover[he_ij.edge()] = 0;
            } else {
                branchCover[he_ij.edge()] = 1;
                total = (total + 1) % 2;
            }*/
        }   
        
        if (singularities[f] != 0 && total == 0) {
            std::cout << "difference at singularity: " << total << std::endl;
        } else if (singularities[f] == 0 && total != 0) {
            std::cout << "difference at non-singularity: " << total << std::endl; 
        }
    }
}

Eigen::SparseMatrix<std::complex<double>> QuadCover::buildLaplacian() {
    size_t n = mesh->nVertices();
    Eigen::SparseMatrix<std::complex<double>> L(n,n);
    std::vector<Eigen::Triplet<std::complex<double>>> triplets;
    
    VertexData<size_t> vertexIndices = mesh->getVertexIndices();
    for (VertexPtr v1 : mesh->vertices()) {
        int index1 = vertexIndices[v1];
        double sum = __DBL_EPSILON__;

        // add neighbor weights
        for (HalfedgePtr heOut : v1.outgoingHalfedges()) {
            VertexPtr v2 = heOut.twin().vertex();
            int index2 = vertexIndices[v2];
            double weight = (geom->cotan(heOut) + geom->cotan(heOut.twin())) / 2;
            sum += weight;

            if (branchCover[heOut] == 1) {
                weight *= -1;
            }
            triplets.push_back(Eigen::Triplet<std::complex<double>>(index1, index2, std::complex<double>(-weight,0)));
        }

        // add diagonal weight
        triplets.push_back(Eigen::Triplet<std::complex<double>>(index1, index1, std::complex<double>(sum,0)));  
    }
    L.setFromTriplets(triplets.begin(), triplets.end());
    return L;
}

VertexData<double> QuadCover::computeOffset() {
    std::cout << "Computing Offset!" << std::endl;
    
    Eigen::SparseMatrix<std::complex<double>> L = buildLaplacian();
    Eigen::SparseMatrix<std::complex<double>> M = assembleM();
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>> solver;
    solver.compute(L);

    // u <- UniformRand(-1,1)
    Eigen::MatrixXcd u = Eigen::MatrixXcd::Random(mesh->nVertices(),1);
    Eigen::MatrixXcd x;

    // inverse power iteration to find eigenvector belonging to the smallest eigenvalue
    for (int i = 0; i < nPowerIterations; i++) {
        x = solver.solve(M * u);
        std::complex<double> norm2 = (x.transpose() * M * x)(0,0);
        u = x / sqrt(norm2);
    }

    // store into VertexData
    VertexData<double> offset(mesh);
    VertexData<size_t> vertexIndices = mesh->getVertexIndices();
    for (VertexPtr v : mesh->vertices()) {
        size_t index = vertexIndices[v];
        offset[v] = x(index,0).real();
    }

    emitTriangles(offset);
    return offset;
}

void QuadCover::emitTriangles(VertexData<double> offsets) {
    std::cout << "Emitting Triangles!" << std::endl;
    
    std::ofstream outfile ("branchcover.obj");
    // write vertices
    for (VertexPtr v : mesh->vertices()) {
        Vector3 pos = geom->position(v) + offsets[v] * geom->normal(v);
        Vector3 neg = geom->position(v) - offsets[v] * geom->normal(v);

        outfile << "v " << pos.x << " " << pos.y << " " << pos.z << std::endl;
        outfile << "v " << neg.x << " " << neg.y << " " << neg.z << std::endl;
    }

    // write face indices
    VertexData<size_t> vertexIndices = mesh->getVertexIndices();
    for (FacePtr f : mesh->faces()) {
        HalfedgePtr he = f.halfedge();

       for (int currSheet = 0; currSheet < 2; currSheet++) {
           outfile << "f ";
           do {
               size_t index = vertexIndices[he.vertex()];
               if (currSheet == 0) {
                   outfile << (2*index)+1 << " ";
               } else {
                   outfile << (2*index+1)+1 << " ";
               }
               currSheet = (currSheet + branchCover[he]) % 2;
               he = he.next();
           } while (he != f.halfedge());
           outfile << std::endl;
       }
    }

    outfile.close();
    std::cout << "Done!" << std::endl;
}
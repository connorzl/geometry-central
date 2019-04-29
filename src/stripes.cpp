#include "geometrycentral/stripes.h"
#include <Eigen/SparseCholesky>

// ONLY WORKS FOR MESHES WITHOUT BOUNDARY
Stripes::Stripes(HalfedgeMesh* m, Geometry<Euclidean>* g) : mesh(m), geom(g), phi(m), r(m), field(m), 
                                                                singularities(m), branchCover(m), omega(m) {
    assert(mesh->nBoundaryLoops() == 0);
}

void Stripes::setup() {
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

Eigen::SparseMatrix<std::complex<double>> Stripes::assembleM() {
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

Eigen::SparseMatrix<std::complex<double>> Stripes::assembleA() {
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

Eigen::MatrixXcd Stripes::principalEigenvector(Eigen::SparseMatrix<std::complex<double>> A, Eigen::SparseMatrix<std::complex<double>> B) {
    // LL^T <- Cholesky(A)
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>> solver;
    solver.compute(A);

    // u <- UniformRand(-1,1)
    Eigen::MatrixXcd x = Eigen::MatrixXcd::Random(mesh->nVertices(),1);

    // inverse power iteration to find eigenvector belonging to the smallest eigenvalue
    for (int i = 0; i < nPowerIterations; i++) {
        x = solver.solve(B * x);
        std::complex<double> norm2 = (x.transpose() * B * x)(0,0);
        x = x / sqrt(norm2);
    }
    return x;
} 

Eigen::MatrixXcd Stripes::principalEigenvector(Eigen::SparseMatrix<double> A, Eigen::SparseMatrix<double> B) {
    // LL^T <- Cholesky(A)
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);

    // u <- UniformRand(-1,1)
    Eigen::MatrixXcd x = Eigen::MatrixXcd::Random(2*mesh->nVertices(),1);

    // inverse power iteration to find eigenvector belonging to the smallest eigenvalue
    for (int i = 0; i < nPowerIterations; i++) {
        x = solver.solve(B * x);
        std::complex<double> norm2 = (x.transpose() * B * x)(0,0);
        x = x / sqrt(norm2);
    }
    return x;
}

VertexData<std::complex<double>> Stripes::computeField() {
    std::cout << "Computing Cross Field... ";
    // Algorithm 1 : Setup
    setup();

    // Algorithm 2 : Smoothest Field
    Eigen::SparseMatrix<std::complex<double>> A = assembleA();
    Eigen::SparseMatrix<std::complex<double>> M = assembleM();
    A = A + eps * M;
    Eigen::MatrixXcd x = principalEigenvector(A,M);

    // map resulting vector to VertexData
    VertexData<size_t> vertexIndices = mesh->getVertexIndices();
    for (VertexPtr v : mesh->vertices()) {
        std::complex<double> c = x(vertexIndices[v],0);
        if (std::abs(c) == 0) {
            field[v] = 0;
        } else {
            field[v] = c / std::abs(c);   
        }
    }
    std::cout << "Done!" << std::endl;
    return field;
}

FaceData<int> Stripes::computeSingularities() {
    std::cout << "Computing Singularities... ";
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
    int total = 0;
    for (FacePtr f : mesh->faces()) {
        double w_ij = w[f.halfedge()];
        double w_jk = w[f.halfedge().next()];
        double w_ki = w[f.halfedge().prev()];
        double Omega_ijk = Omega[f];
        double phi = (w_ij + w_jk + w_ki - Omega_ijk) / (2.0 * M_PI);
        singularities[f] = std::round(phi);
        total += singularities[f];
    }
    std::cout << "Sum: " << total << std::endl;
    return singularities;
}

void Stripes::edgeData() {
    std::cout<< "Computing EdgeData... ";
    std::complex<double> i(0, 1);
    for (EdgePtr e : mesh->edges()) {
        HalfedgePtr he_ij = e.halfedge();
        HalfedgePtr he_ji = he_ij.twin();

        // disambiguate the big vectors
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
        
        double sign;
        if ( ang >= -M_PI_2 && ang < M_PI_2 ) {
            sign = 1;
            branchCover[e] = 0;
        } else {
            sign = -1;
            branchCover[e] = 1;
        }
    
        double phi_i = std::arg(f_ij);
        double phi_j = std::arg(sign * f_ji);
        double l_ij = geom->length(e);
        omega[e] = lambda * (l_ij / 2.0) * (cos(phi_i - theta_ij) + cos(phi_j - theta_ji));
    }
    std::cout << "Done!" << std::endl;
}

Eigen::SparseMatrix<double> Stripes::EnergyMatrix() {
    size_t n = mesh->nVertices();
    Eigen::SparseMatrix<double> A(2*n,2*n);
    std::vector<Eigen::Triplet<double>> triplets;

    VertexData<size_t> vertexIndices = mesh->getVertexIndices();
    for (EdgePtr e : mesh->edges()) {
        HalfedgePtr he_ij = e.halfedge();
        
        // cotan weights
        double cotA = geom->cotan(he_ij);
        double cotB = geom->cotan(he_ij.twin());
        if (singularities[he_ij.face()] != 0) cotA = 0;
        if (singularities[he_ij.twin().face()] != 0) cotB = 0;
        double w = (cotA + cotB) / 2.0;

        // indices
        size_t i = 2 * vertexIndices[he_ij.vertex()];
        size_t j = 2 * vertexIndices[he_ij.twin().vertex()];

        // add diagonal terms
        triplets.push_back( Eigen::Triplet<double>(i,i,w) );
        triplets.push_back( Eigen::Triplet<double>(i+1,i+1,w) );
        triplets.push_back( Eigen::Triplet<double>(j,j,w) );
        triplets.push_back( Eigen::Triplet<double>(j+1,j+1,w) );

        // transport coefficient components
        double x = w * cos(omega[e]);
        double y = w * sin(omega[e]);

        // stays on same sheet
        if (branchCover[e] == 0) {
            // A_ij
            triplets.push_back( Eigen::Triplet<double>(i,j,-x) ); triplets.push_back( Eigen::Triplet<double>(i,j+1,-y) );
            triplets.push_back( Eigen::Triplet<double>(i+1,j,y) ); triplets.push_back( Eigen::Triplet<double>(i+1,j+1,-x) );
            // A_ji
            triplets.push_back( Eigen::Triplet<double>(j,i,-x) ); triplets.push_back( Eigen::Triplet<double>(j,i+1,y) );
            triplets.push_back( Eigen::Triplet<double>(j+1,i,-y) ); triplets.push_back( Eigen::Triplet<double>(j+1,i+1,-x) );
        } else {
            // A_ij
            triplets.push_back( Eigen::Triplet<double>(i,j,-x) ); triplets.push_back( Eigen::Triplet<double>(i,j+1,y) );
            triplets.push_back( Eigen::Triplet<double>(i+1,j,y) ); triplets.push_back( Eigen::Triplet<double>(i+1,j+1,x) );
            // A_ji
            triplets.push_back( Eigen::Triplet<double>(j,i,-x) ); triplets.push_back( Eigen::Triplet<double>(j,i+1,y) );
            triplets.push_back( Eigen::Triplet<double>(j+1,i,y) ); triplets.push_back( Eigen::Triplet<double>(j+1,i+1,x) );
        }
    }

    A.setFromTriplets(triplets.begin(),triplets.end());
    return A;
}

Eigen::SparseMatrix<double> Stripes::MassMatrix() {
    size_t n = mesh->nVertices();
    Eigen::SparseMatrix<double> M(2*n,2*n);
    std::vector<Eigen::Triplet<double>> triplets;

    VertexData<size_t> vertexIndices = mesh->getVertexIndices();
    for (FacePtr f : mesh->faces()) {
        HalfedgePtr he_ij = f.halfedge();
        size_t i = vertexIndices[he_ij.vertex()];
        size_t j = vertexIndices[he_ij.next().vertex()];
        size_t k = vertexIndices[he_ij.prev().vertex()];

        double area = geom->area(f);
        triplets.push_back(Eigen::Triplet<double>(2*i, 2*i, area/3.));
        triplets.push_back(Eigen::Triplet<double>(2*i+1, 2*i+1, area/3.));
        triplets.push_back(Eigen::Triplet<double>(2*j, 2*j, area/3.));
        triplets.push_back(Eigen::Triplet<double>(2*j+1, 2*j+1, area/3.));
        triplets.push_back(Eigen::Triplet<double>(2*k, 2*k, area/3.));
        triplets.push_back(Eigen::Triplet<double>(2*k+1, 2*k+1, area/3.));
    }
    M.setFromTriplets(triplets.begin(),triplets.end());
    return M;
}

void Stripes::computeStripes() {
    std::cout << "Computing Stripes... ";
    Eigen::SparseMatrix<double> A = EnergyMatrix();
    Eigen::SparseMatrix<double> B = MassMatrix();
    Eigen::MatrixXcd x = principalEigenvector(A,B);
    std::cout << "Done!" << std::endl;
    std::cout << x.rows() << "," << x.cols() << std::endl;
}
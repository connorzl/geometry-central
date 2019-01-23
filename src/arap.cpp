#include "geometrycentral/arap.h"
#include<Eigen/SparseCholesky>
#include <Eigen/SVD>
#include <Eigen/Dense>
using namespace geometrycentral;

ARAP::ARAP(HalfedgeMesh* m, Geometry<Euclidean>* g) : mesh(m), geom(g), vertexIndices(mesh), isoTriangleParam(mesh), uvCoords(mesh) {
    vertexIndices = mesh->getVertexIndices();
}

Eigen::SparseMatrix<std::complex<double>> ARAP::createLaplaceMatrix() {
    size_t n = mesh->nVertices();
    Eigen::SparseMatrix<std::complex<double>> A(n,n);
    std::vector<Eigen::Triplet<std::complex<double>>> triplets;
    
    for (VertexPtr v1 : mesh->vertices()) {
        int index1 = vertexIndices[v1];
        double sum = __DBL_EPSILON__;

        // add neighbor weights
        for (HalfedgePtr heOut : v1.outgoingHalfedges()) {
            VertexPtr v2 = heOut.twin().vertex();
            int index2 = vertexIndices[v2];
            double weight = (geom->cotan(heOut) + geom->cotan(heOut.twin())) / 2.0;
            
            sum += weight;
            triplets.push_back(Eigen::Triplet<std::complex<double>>(index1, index2, std::complex<double>(-weight,0)));
        }

        // add diagonal weight
        triplets.push_back(Eigen::Triplet<std::complex<double>>(index1, index1, std::complex<double>(sum,0)));  
    }

    A.setFromTriplets(triplets.begin(), triplets.end());
    return A;
}

void ARAP::computeIsoTriangleParam() {
    for (FacePtr f : mesh->faces()) {
        // Gather elements
        HalfedgePtr he = f.halfedge();
        double l_ab = geom->length(he.edge());
        double l_ac = geom->length(he.prev().edge());
        double theta_a = geom->angle(he.next()); // radians

        // Place first vertex at (0,0)
        isoTriangleParam[he] = Vector2{0,0};

        // Place second vertex at (|ab|,0)
        isoTriangleParam[he.next()] = Vector2{l_ab,0};

        // Place third vertex at (|ac|,0) rotated by theta_a CCW
        isoTriangleParam[he.prev()] = Vector2{cos(theta_a) * l_ac, sin(theta_a) * l_ac}; 
    }
}

FaceData<Eigen::Matrix2d> ARAP::computeLRotations(VertexData<Vector2> const &u) {
    FaceData<Eigen::Matrix2d> L(mesh);
    for (FacePtr f : mesh->faces()) {
        // Gather elements
        HalfedgePtr he1 = f.halfedge();
        HalfedgePtr he2 = he1.next();
        HalfedgePtr he3 = he1.prev();
        std::vector<Vector2> ut = { u[he1.vertex()], u[he2.vertex()], u[he3.vertex()] };
        std::vector<Vector2> xt = { isoTriangleParam[he1], isoTriangleParam[he2], isoTriangleParam[he3] };
        std::vector<double> thetat = { geom->cotan(he1), geom->cotan(he2), geom->cotan(he3) };

        // Compute St matrix
        Eigen::Matrix2d St = Eigen::Matrix2d::Zero();
        for (int i = 0; i < 3; i++) {
            Vector2 ui = ut[i] - ut[(i+1) % 3];
            Vector2 xi = xt[i] - xt[(i+1) % 3];
            
            St(0,0) += thetat[i] * ui.x * xi.x;
            St(0,1) += thetat[i] * ui.x * xi.y;
            St(1,0) += thetat[i] * ui.y * xi.x;
            St(1,1) += thetat[i] * ui.y * xi.y;
        }

        // Perform SVD decomposition, where L_t = UV^T
        Eigen::JacobiSVD<Eigen::Matrix2d> svd( St, Eigen::ComputeFullU | Eigen::ComputeFullV );
        Eigen::Matrix2d U = svd.matrixU();
        Eigen::Matrix2d V = svd.matrixV();

        Eigen::Matrix2d UVT = U * V.transpose();
        if (UVT.determinant() < 0) {
            V.col(1) *= -1;
            UVT = U * V.transpose();
        }
        L[f] = UVT;
    }
    return L;
}

Eigen::MatrixXcd ARAP::computebVector(FaceData<Eigen::Matrix2d> const &L) {
    Eigen::MatrixXcd b = Eigen::MatrixXcd::Zero(mesh->nVertices(),1);
    for (VertexPtr v : mesh->vertices()) {
        size_t index = vertexIndices[v];

        for (HalfedgePtr he_ij : v.outgoingHalfedges()) {
            HalfedgePtr he_ji = he_ij.twin();

            // first triangle term
            if (he_ij.isReal()) {
                Vector2 xi = isoTriangleParam[he_ij];
                Vector2 xj = isoTriangleParam[he_ij.next()];

                double cotan_ij = geom->cotan(he_ij);
                Eigen::Matrix2d Lt_ij = L[he_ij.face()];

                std::complex<double> sub((xi-xj).x, (xi-xj).y); 
                std::complex<double> rot(Lt_ij(0,0), Lt_ij(1,0));
                b(index,0) += cotan_ij * rot * sub / 2.0;
            }

            // second triangle term
            if (he_ji.isReal()) {
                Vector2 xi = isoTriangleParam[he_ji.next()];
                Vector2 xj = isoTriangleParam[he_ji];

                double cotan_ji = geom->cotan(he_ji);
                Eigen::Matrix2d Lt_ji = L[he_ji.face()];  

                std::complex<double> sub((xi-xj).x, (xi-xj).y);    
                std::complex<double> rot(Lt_ji(0,0), Lt_ji(1,0));
                b(index,0) += cotan_ji * rot * sub / 2.0;
            }
        } 
    }
    return b;
}

void ARAP::computeARAP() {
    // Build Laplace Matrix A (n x n) and factorize
    Eigen::SparseMatrix<std::complex<double>> A = createLaplaceMatrix();
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>> solver;
    solver.compute(A);

    // Compute isometric parameterization for each triangle t
    computeIsoTriangleParam();

    // Initial parameterization u (using SCP)
    SpectralConformal s = SpectralConformal(mesh,geom);
    VertexData<Vector2> u = s.computeSpectralConformal();
    
    // Repeat the following until convergence:
    for (int i = 0; i < 10; i++) {
        // Fix the mapping u (n x 1) and solve for L_t (2x2) for each triangle t
        FaceData<Eigen::Matrix2d> L = computeLRotations(u);

        // Compute b (n x 1) using L
        Eigen::MatrixXcd b = computebVector(L);

        // Solve Au = b
        Eigen::MatrixXcd u_new = solver.solve(b);

        // Update u
        for (VertexPtr v : mesh->vertices()) {
            std::complex<double> uv = u_new(vertexIndices[v],0);
            u[v] = Vector2{uv.real(), uv.imag()};
        }
        std::cout << "finished iteration: " << i << std::endl;
    }

    // normalize
    uvCoords = u;
    normalize();

    // write output obj file
    std::ofstream outfile ("ARAP.obj");
    writeToFile(outfile);
    outfile.close();
    std::cout<<"Done ARAP!"<<std::endl;
}

void ARAP::writeToFile(std::ofstream &outfile) {
    // write vertices
    for (VertexPtr v : mesh->vertices()) {
        outfile << "v " << geom->position(v).x << " " << geom->position(v).y << " " << geom->position(v).z << std::endl;
    }

    // write uvs
    for (VertexPtr v : mesh->vertices()) {
        outfile << "vt " << uvCoords[v].x << " " << uvCoords[v].y << std::endl;
    }

    // write indices
    VertexData<size_t> index = mesh->getVertexIndices();
    for (FacePtr f : mesh->faces()) {
       HalfedgePtr he = f.halfedge();
       outfile << "f";
       do {
           VertexPtr v = he.vertex();
           outfile << " " << index[v] + 1 << "/" << index[v] + 1;

           he = he.next();
       } while (he != f.halfedge());
       outfile << std::endl;
    }

   outfile.close();
}

void ARAP::normalize() {
    // compute center of mass
    Vector2 cm = {0,0};
    for (VertexPtr v : mesh->vertices()) {
        Vector2 uv = uvCoords[v];
        cm += uv;
    }
    cm /= mesh->nVertices();

    double r = 0;
    for (VertexPtr v : mesh->vertices()) {
        Vector2 &uv = uvCoords[v];
        uv -= cm;
        r = std::max(r, norm(uv));
    }

    for (VertexPtr v : mesh->vertices()) {
        Vector2 &uv = uvCoords[v];
        uv /= r;
    }
}
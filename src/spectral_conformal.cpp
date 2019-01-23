#include "geometrycentral/spectral_conformal.h"
#include<Eigen/SparseCholesky>

using namespace geometrycentral;

SpectralConformal::SpectralConformal(HalfedgeMesh* m, Geometry<Euclidean>* g) : mesh(m), geom(g), vertexIndices(mesh), uvCoords(mesh) {}

Eigen::SparseMatrix<std::complex<double>> SpectralConformal::createEDMatrix() {
    size_t n = mesh->nVertices();
    Eigen::SparseMatrix<std::complex<double>> ED(n,n);
    std::vector<Eigen::Triplet<std::complex<double>>> triplets;
    
    for (VertexPtr v1 : mesh->vertices()) {
        int index1 = vertexIndices[v1];
        double sum = __DBL_EPSILON__;

        // add neighbor weights
        for (HalfedgePtr heOut : v1.outgoingHalfedges()) {
            VertexPtr v2 = heOut.twin().vertex();
            int index2 = vertexIndices[v2];
            double weight = (geom->cotan(heOut) + geom->cotan(heOut.twin())) / 2;
            
            sum += weight;
            triplets.push_back(Eigen::Triplet<std::complex<double>>(index1, index2, std::complex<double>(-weight,0)));
        }

        // add diagonal weight
        triplets.push_back(Eigen::Triplet<std::complex<double>>(index1, index1, std::complex<double>(sum,0)));  
    }

    ED.setFromTriplets(triplets.begin(), triplets.end());
    return ED;
}

Eigen::SparseMatrix<std::complex<double>> SpectralConformal::createAMatrix() {
    size_t n = mesh->nVertices();
    Eigen::SparseMatrix<std::complex<double>> A(n,n);
    std::complex<double> c1 (0,0.25);
    std::complex<double> c2 (0,-0.25);
    std::vector<Eigen::Triplet<std::complex<double>>> triplets;
    for (HalfedgePtr he : mesh->imaginaryHalfedges()) {
        size_t i = vertexIndices[he.vertex()];
        size_t j = vertexIndices[he.twin().vertex()];

        triplets.push_back(Eigen::Triplet<std::complex<double>>(i, j, c1));
        triplets.push_back(Eigen::Triplet<std::complex<double>>(j, i, c2));
    }
    A.setFromTriplets(triplets.begin(),triplets.end());
    return A;
}

// Compute Residual = Ax − (x^T Ax)x, where A is n x n and x is n x 1 unit vector
double SpectralConformal::computeResidual(Eigen::SparseMatrix<std::complex<double>> A, Eigen::MatrixXcd x) {
    Eigen::MatrixXcd Ax = A * x;                       // n x 1
    Eigen::MatrixXcd xH = x.adjoint();                 // 1 x n
    Eigen::MatrixXcd xHAx = xH * Ax;                   // 1 x 1
    Eigen::MatrixXcd xHx = xH * x;                     // 1 x 1
    std::complex<double> lambda = xHAx(0,0) / xHx(0,0);

    Eigen::MatrixXcd resid = Ax - lambda * x;          // n x 1
    double r = resid.norm();

    std::cout<< r <<std::endl;
    return r;
}

VertexData<Vector2> SpectralConformal::computeSpectralConformal() {
    // compute vertex indices
    vertexIndices = mesh->getVertexIndices();
    
    // Build EC matrix = Dirichlet energy matrix - area matrix
    Eigen::SparseMatrix<std::complex<double>> ED = 0.5 * createEDMatrix();
    Eigen::SparseMatrix<std::complex<double>> A = createAMatrix();
    Eigen::SparseMatrix<std::complex<double>> EC = ED - A;

    // Solve using inverse power method 
    size_t n = EC.rows();
    Eigen::SimplicialLLT<Eigen::SparseMatrix<std::complex<double>>> solver;
    solver.compute(EC);
    Eigen::MatrixXcd ones = Eigen::MatrixXcd::Ones(n,1);
    Eigen::MatrixXcd x = Eigen::MatrixXcd::Random(n,1);

    do {
        // solve step
        x = solver.solve(x);

        // subtract mean
        std::complex<double> mean = x.sum() / (double)n;
        x = x - mean * ones;
        
        // normalize
        x = x / x.norm();

    } while (computeResidual(EC,x) > 1.0f * pow(10,-10));

    // assign flattening
    for (VertexPtr v : mesh->vertices()) {
        size_t index = vertexIndices[v];
        std::complex<double> c = x(index,0);
        uvCoords[v] = Vector2{c.real(), c.imag()};
    }

    // normalize
    normalize();

    // write output obj file
    std::ofstream outfile ("SpectralConformal.obj");
    writeToFile(outfile);
    outfile.close();
    std::cout<<"Done SCP!"<<std::endl;

    return uvCoords;
}

void SpectralConformal::writeToFile(std::ofstream &outfile) {
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

void SpectralConformal::normalize() {
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
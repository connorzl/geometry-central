#include "geometrycentral/tutte.h"
#include<Eigen/SparseCholesky>

using namespace geometrycentral;

Tutte::Tutte(HalfedgeMesh* m, Geometry<Euclidean>* g) : mesh(m), geom(g), vertexIndices(mesh), uvCoords(mesh) {}

void Tutte::separateVertices() {
     // remap the vertex indices to interior and boundary
    int iN = 0;
    int bN = mesh->nInteriorVertices();

    for (VertexPtr v : mesh->vertices()) {
        if (!v.isBoundary()) {
            vertexIndices[v] = iN++;
        } else {
            vertexIndices[v] = bN++;
        }
    }
}

Eigen::MatrixXd Tutte::mapToUnitCircle() {
    int iN = mesh->nInteriorVertices();
    int bN = mesh->nImaginaryHalfedges();
    Eigen::MatrixXd P1(bN, 2);

    int i = 0;
    float doublePi = 2.0f * M_PI;

    for (HalfedgePtr he : mesh->imaginaryHalfedges()) {
        VertexPtr v = he.vertex();
        int index = vertexIndices[v] - iN;

        P1(index, 0) = cos( doublePi * i / bN );
        P1(index, 1) = sin( doublePi * i / bN );
        i++;
    }

    return P1;
}

Eigen::SparseMatrix<double> Tutte::createKMatrix() {
    size_t n = mesh->nVertices();
    Eigen::SparseMatrix<double> K(n,n);
    std::vector<Eigen::Triplet<double>> triplets;
    
    for (VertexPtr v1 : mesh->vertices()) {
        int index1 = vertexIndices[v1];
        double sum = __DBL_EPSILON__;

        // add neighbor weights
        for (HalfedgePtr heOut : v1.outgoingHalfedges()) {
            VertexPtr v2 = heOut.twin().vertex();
            int index2 = vertexIndices[v2];
            double weight = std::fmax(0,(geom->cotan(heOut) + geom->cotan(heOut.twin())) / 2);
            
            sum += weight;
            triplets.push_back(Eigen::Triplet<double>(index1, index2, -weight));    
        }

        // add diagonal weight
        triplets.push_back(Eigen::Triplet<double>(index1, index1,sum));  
    }

    K.setFromTriplets(triplets.begin(), triplets.end());
    return K;
}

void Tutte::writeToFile(std::ofstream &outfile) {
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

void Tutte::normalize() {
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

void Tutte::computeTutteEmbedding() {
    // compute indices where interior vertices come before boundary vertices
    separateVertices();

    // Map exterior vertices to unit circle (P1)
    Eigen::MatrixXd P2 = mapToUnitCircle();

    // Build K matrix
    Eigen::SparseMatrix<double> K = createKMatrix();

    // Extract K1, B^T, B, K2
    int iN = mesh->nInteriorVertices();
    int bN = mesh->nImaginaryHalfedges();
    Eigen::SparseMatrix<double> K1 = K.block(0, 0, iN, iN); 
    Eigen::SparseMatrix<double> BT = K.block(0, iN, iN, bN);

    // Solve for interior vertex uv-coordinates (P2)
    // K1 * P1 + BT * P2 = 0 -> P1 = (K1)^-1 * -BT * P2 
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(K1);
    if(solver.info()!= Eigen::Success) {
        std::cout << "Decomposition Failed!" << std::endl;
        return;
    }   

    Eigen::MatrixXd P1 = solver.solve(-BT * P2);
    if(solver.info()!= Eigen::Success) {
        std::cout << "Solving Failed!" << std::endl;
        return;
    }   

    // store uv coordinates
    for (VertexPtr v: mesh->vertices()) {
        int index = vertexIndices[v];
        if (v.isBoundary()) {
            index -= iN;
            uvCoords[v] = Vector2{P2(index, 0), P2(index, 1)};
        } else {
            uvCoords[v] = Vector2{P1(index, 0), P1(index, 1)};
        }
    }

    // normalize
    normalize();

    // write output obj file
    std::ofstream outfile ("tutte.obj");
    writeToFile(outfile);
    outfile.close();
    std::cout<<"Done!"<<std::endl;
}
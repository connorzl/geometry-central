#include "geometrycentral/tutte.h"
#include<Eigen/SparseCholesky>

using namespace geometrycentral;

Tutte::Tutte(HalfedgeMesh* m, Geometry<Euclidean>* g) : mesh(m), geom(g), vertexIndices(mesh), uvCoords(mesh) {}

void Tutte::separateVertices() {
     // separate the vertices of the longest boundary loop from other boundary vertices
    /*
    size_t longestIndex = mesh->longestBoundaryLoop();
    for (size_t i = 0; i < mesh->nBoundaryLoops(); i++) {
        BoundaryPtr b = mesh->boundaryLoop(i);
        HalfedgePtr he = b.halfedge();
        HalfedgePtr curr = he;
        do {
            assert(!curr.isReal() && curr.vertex().isBoundary());
            
            if (i == longestIndex) {
                boundaryVertices.push_back(curr.vertex());
            } else {
                interiorVertices.push_back(curr.vertex());
            }
            curr = curr.next();
        } while(curr != he);
    }
    */
    // now add in the interior vertices
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

    for (HalfedgePtr he: mesh->imaginaryHalfedges()) {
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
    
/*
    for (size_t i = 0; i < n; i++) {
        VertexPtr v1;
        if ( i < boundaryVertices.size() ) {
            v1 = boundaryVertices[i];
        } else {
            v1 = interiorVertices[i - boundaryVertices.size()];
        }
 */   
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
    // Extract interior and exterior vertices
    separateVertices();

    // Map exterior vertices to unit circle (P1)
    Eigen::MatrixXd P1 = mapToUnitCircle();

    // Build K matrix
    Eigen::SparseMatrix<double> K = createKMatrix();

    
    // Extract K1, B^T, B, K2
    int iN = mesh->nInteriorVertices();
    int bN = mesh->nImaginaryHalfedges();
    Eigen::SparseMatrix<double> B = K.block(0, 0, iN, iN); // r: 0 - iN, 0 - iN
    Eigen::SparseMatrix<double> K2 = K.block(0, iN, iN, bN); // r: 0 - iN, c: iN - N


    // Solve for interior vertex uv-coordinates (P2)
    // B * P1 + K2 * P2 = 0 -> P2 = (K2)^-1 * -B * P1 -> flipped
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver;
    solver.compute(B);
    if(solver.info()!= Eigen::Success) {
        std::cout << "Decomposition Failed!" << std::endl;
        return;
    }   

    Eigen::MatrixXd P2 = solver.solve(-K2 * P1);
    if(solver.info()!= Eigen::Success) {
        std::cout << "Solving Failed!" << std::endl;
        return;
    }   

    // store uv coordinates
    for (VertexPtr v: mesh->vertices()) {
        int index = vertexIndices[v];

        if (v.isBoundary()) {
            index -= iN;
            uvCoords[v] = Vector2{P1(index, 0), P1(index, 1)};

        } else {
            uvCoords[v] = Vector2{P2(index, 0), P2(index, 1)};
        }
    }

        std::cout<<"got this far 1"<<std::endl;

    // normalize
    normalize();

    // write output obj file
    std::ofstream outfile ("tutte.obj");
    writeToFile(outfile);
    outfile.close();

    // sanity check for storage
    //assert(mesh->getVertexIndices()[interiorVertices[0]] == 0);
    //assert(P2(0,0) == geom->uvCoords[mesh->vertex(0)].x);
    //assert(P2(0,1) == geom->uvCoords[mesh->vertex(0)].y);
    std::cout<<"Done!"<<std::endl;
}
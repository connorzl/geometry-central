#include "geometrycentral/least_squares_conformal.h"
#include<Eigen/SparseCholesky>

using namespace geometrycentral;

LSCM::LSCM(HalfedgeMesh* m, Geometry<Euclidean>* g) : mesh(m), geom(g), vertexIndices(mesh), uvCoords(mesh) {}

void LSCM::separateVertices(VertexPtr v1, VertexPtr v2) {
     // remap the vertex indices to interior and boundary
    int iN = 0;
    int bN = mesh->nVertices()-2;

    for (VertexPtr v : mesh->vertices()) {
        if (v == v1 || v == v2) {
            vertexIndices[v] = bN++;
        } else {
            vertexIndices[v] = iN++;
        }
    }
}

Eigen::SparseMatrix<std::complex<double>> LSCM::createEDMatrix() {
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

Eigen::SparseMatrix<std::complex<double>> LSCM::createAMatrix() {
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

VertexData<Vector2> LSCM::computeLSCM() {
    // Find 2 fixed vertices on the boundary
    VertexPtr v1;
    VertexPtr v2;
    double maxDistance = 0;
    for (HalfedgePtr p1 : mesh->imaginaryHalfedges()) {
        for (HalfedgePtr p2 : mesh->imaginaryHalfedges()) {
            double distance = norm(geom->position(p1.vertex()) - geom->position(p2.vertex()));
            if (distance > maxDistance) {
                v1 = p1.vertex();
                v2 = p2.vertex();
                maxDistance = distance;
            }
        }
    }
    // fix these two vertices
    Eigen::MatrixXcd P2(2,1);
    P2(0,0) = std::complex<double>(0,-1);
    P2(1,0) = std::complex<double>(0,1);
    
    // compute vertex indices, separate vertices into [ unknown, fixed ]^T
    separateVertices(v1,v2);
    
    // Build EC matrix = Dirichlet energy matrix - area matrix
    Eigen::SparseMatrix<std::complex<double>> ED = 0.5 * createEDMatrix();
    Eigen::SparseMatrix<std::complex<double>> A = createAMatrix();
    Eigen::SparseMatrix<std::complex<double>> EC = ED - A;

    // Extract K1, B^T, B, K2
    int iN = mesh->nVertices()-2;
    int bN = 2;
    Eigen::SparseMatrix<std::complex<double>> K1 = EC.block(0, 0, iN, iN); 
    Eigen::SparseMatrix<std::complex<double>> BT = EC.block(0, iN, iN, bN);

    // Solve for interior vertex uv-coordinates (P1)
    // K1 * P1 + BT * P2 = 0 -> K1 * P1 = -BT * P2 
    Eigen::SimplicialLLT<Eigen::SparseMatrix<std::complex<double>>> solver;
    solver.compute(K1);
    if(solver.info()!= Eigen::Success) {
        std::cout << "Decomposition Failed!" << std::endl;
        return nullptr;
    }   

    Eigen::MatrixXcd P1 = solver.solve(-BT * P2);
    if(solver.info()!= Eigen::Success) {
        std::cout << "Solving Failed!" << std::endl;
        return nullptr;
    }   

    // assign flattening
    for (VertexPtr v : mesh->vertices()) {
        int index = vertexIndices[v];
        if (v == v1 || v == v2) {
            index -= iN;
            uvCoords[v] = Vector2{P2(index, 0).real(), P2(index, 0).imag()};
        } else {
            uvCoords[v] = Vector2{P1(index, 0).real(), P1(index, 0).imag()};
        }
    }

    // normalize
    normalize();

    // write output obj file
    std::ofstream outfile ("LSCM.obj");
    writeToFile(outfile);
    outfile.close();
    std::cout<<"Done LSCM!"<<std::endl;

    return uvCoords;
}

void LSCM::writeToFile(std::ofstream &outfile) {
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

void LSCM::normalize() {
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
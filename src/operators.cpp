#include "geometrycentral/operators.h"

Eigen::SparseMatrix<double> Operators::laplaceMatrix(HalfedgeMesh* mesh, Geometry<Euclidean>* geom) {
    size_t n = mesh->nVertices();
    Eigen::SparseMatrix<double> L(n,n);
    std::vector<Eigen::Triplet<double>> triplets;
    VertexData<size_t> vertexIndices = mesh->getVertexIndices();
    
    for (VertexPtr v1 : mesh->vertices()) {
        int index1 = vertexIndices[v1];
        double sum = 1e-8;

        // add neighbor weights
        for (HalfedgePtr heOut : v1.outgoingHalfedges()) {
            VertexPtr v2 = heOut.twin().vertex();
            int index2 = vertexIndices[v2];
            double weight = (geom->cotan(heOut) + geom->cotan(heOut.twin())) / 2;
            
            sum += weight;
            triplets.push_back(Eigen::Triplet<double>(index1, index2,-weight));
        }

        // add diagonal weight
        triplets.push_back(Eigen::Triplet<double>(index1,index1,sum));  
    }

    L.setFromTriplets(triplets.begin(), triplets.end());
    return L;
}

HalfedgeData<double> Operators::computeAngles(HalfedgeMesh* mesh, EdgeData<double> &edgeLengths) {
    HalfedgeData<double> angles(mesh);
    for (FacePtr f : mesh->faces()) {
        HalfedgePtr he_ij = f.halfedge();
        HalfedgePtr he_jk = he_ij.next();
        HalfedgePtr he_ki = he_ij.prev();

        double l_ij = edgeLengths[he_ij.edge()];
        double l_jk = edgeLengths[he_jk.edge()];
        double l_ki = edgeLengths[he_ki.edge()];

        angles[he_ij] = acos( (pow(l_ij,2) - pow(l_jk,2) - pow(l_ki,2)) / (-2 * l_jk * l_ki) );
        angles[he_jk] = acos( (pow(l_jk,2) - pow(l_ij,2) - pow(l_ki,2)) / (-2 * l_ij * l_ki) );
        angles[he_ki] = acos( (pow(l_ki,2) - pow(l_ij,2) - pow(l_jk,2)) / (-2 * l_ij * l_jk) );

        if (std::isnan(angles[he_ij]) || std::isnan(angles[he_jk]) || std::isnan(angles[he_ki])) {
            std::cout << (pow(l_ij,2) - pow(l_jk,2) - pow(l_ki,2)) / (-2 * l_jk * l_ki) << std::endl;
            std::cout << (pow(l_jk,2) - pow(l_ij,2) - pow(l_ki,2)) / (-2 * l_ij * l_ki) << std::endl;
            std::cout << (pow(l_ki,2) - pow(l_ij,2) - pow(l_jk,2)) / (-2 * l_ij * l_jk) << std::endl;
            std::cout << l_ij << "," << l_jk << "," << l_ki << std::endl;
            throw std::logic_error ("triangle inequality violated when computing angles");
        }
    }
    return angles;
}

Eigen::SparseMatrix<double> Operators::intrinsicLaplaceMatrix(HalfedgeMesh* mesh, EdgeData<double> &edgeLengths) {
    HalfedgeData<double> angles = Operators::computeAngles(mesh, edgeLengths);

    size_t n = mesh->nVertices();
    Eigen::SparseMatrix<double> L(n,n);
    std::vector<Eigen::Triplet<double>> triplets;
    VertexData<size_t> vertexIndices = mesh->getVertexIndices();
    for (VertexPtr v1 : mesh->vertices()) {
        int index1 = vertexIndices[v1];
        double sum = 0;

        // add neighbor weights
        for (HalfedgePtr heOut : v1.outgoingHalfedges()) {
            VertexPtr v2 = heOut.twin().vertex();
            int index2 = vertexIndices[v2];
            double cot_ij, cot_ji, weight;
            if (heOut.isReal()) {
                cot_ij = 1.0 / tan(angles[heOut]);
            } else {
                cot_ij = 0;
            }
            if (heOut.twin().isReal()) {
                cot_ji = 1.0 / tan(angles[heOut.twin()]);
            } else {
                cot_ji = 0;
            }
            weight = (cot_ij + cot_ji) / 2; 
            
            sum += weight;
            triplets.push_back(Eigen::Triplet<double>(index1, index2,-weight));
        }

        // add diagonal weight
        triplets.push_back(Eigen::Triplet<double>(index1,index1,sum));  
    }
    L.setFromTriplets(triplets.begin(), triplets.end());
    return L;
}

Eigen::MatrixXd Operators::intrinsicCurvature(HalfedgeMesh* mesh, EdgeData<double> &edgeLengths) {
    HalfedgeData<double> theta = Operators::computeAngles(mesh, edgeLengths);
    
    size_t n = mesh->nVertices();
    Eigen::MatrixXd K(n,1);
    VertexData<size_t> vertexIndices = mesh->getVertexIndices();
    for (VertexPtr v : mesh->vertices()) {
        double sum = 0;
        for (HalfedgePtr he : v.outgoingHalfedges()) {
            if (he.edge().isBoundary()) continue;
            sum += theta[he.next()];
        }
        if (v.isBoundary()) {
            K(vertexIndices[v],0) = M_PI - sum; // this doesn't really matter since u is set to 0
        } else {
            K(vertexIndices[v],0) = 2*M_PI - sum;
        }
    }
    return K;
}

Eigen::SparseMatrix<double> Operators::buildExteriorDerivative0Form(HalfedgeMesh* mesh, Geometry<Euclidean>* geom) {
    Eigen::SparseMatrix<double> d0(mesh->nEdges(),mesh->nVertices());
    EdgeData<size_t> edgeIndices = mesh->getEdgeIndices();
    VertexData<size_t> vertexIndices = mesh->getVertexIndices();
    std::vector<Eigen::Triplet<double>> triplets;
    for (EdgePtr e : mesh->edges()) {
	    VertexPtr v1 = e.halfedge().vertex();
        VertexPtr v2 = e.halfedge().twin().vertex();
        triplets.push_back(Eigen::Triplet<double>(edgeIndices[e],vertexIndices[v2],1));
        triplets.push_back(Eigen::Triplet<double>(edgeIndices[e],vertexIndices[v1],-1));
    }
    d0.setFromTriplets(triplets.begin(), triplets.end());
    return d0;
}

Eigen::SparseMatrix<double> Operators::buildHodgeStar1Form(HalfedgeMesh* mesh, Geometry<Euclidean>* geom) {
    Eigen::SparseMatrix<double> hodge1(mesh->nEdges(),mesh->nEdges());
    std::vector<Eigen::Triplet<double>> triplets;
    EdgeData<size_t> edgeIndices = mesh->getEdgeIndices();
    for (EdgePtr e : mesh->edges()) {
        double cotA = geom->cotan(e.halfedge());
        double cotB = geom->cotan(e.halfedge().twin());
        double ratio = 0.5 * (cotA + cotB);
        triplets.push_back(Eigen::Triplet<double>(edgeIndices[e],edgeIndices[e],ratio));
    }
    hodge1.setFromTriplets(triplets.begin(), triplets.end());
    return hodge1;
}

void Operators::hyperbolicEdgeFlips(HalfedgeMesh* mesh, EdgeData<double> &edgeLengthsCM) {
    int numFlips;
    do {
        numFlips = 0;
        for (EdgePtr e : mesh->edges()) {
            if (e.isBoundary()) continue;
            HalfedgePtr he_ij = e.halfedge();
            HalfedgePtr he_ji = he_ij.twin();

            double l_alpha_ij = edgeLengthsCM[he_ij.edge()];
            double l_beta_ij = edgeLengthsCM[he_ij.next().edge()];
            double l_gamma_ij = edgeLengthsCM[he_ij.prev().edge()];

            double l_alpha_ji = edgeLengthsCM[he_ji.edge()];
            double l_beta_ji = edgeLengthsCM[he_ji.prev().edge()];
            double l_gamma_ji = edgeLengthsCM[he_ji.next().edge()];

            double alpha_ij = l_alpha_ij / (l_beta_ij * l_gamma_ij);
            double beta_ij = l_beta_ij / (l_alpha_ij * l_gamma_ij);
            double gamma_ij = l_gamma_ij / (l_alpha_ij * l_beta_ij);

            double alpha_ji = l_alpha_ji / (l_beta_ji * l_gamma_ji);
            double beta_ji = l_beta_ji / (l_alpha_ji * l_gamma_ji);
            double gamma_ji = l_gamma_ji / (l_alpha_ji * l_beta_ji);

            double sum = beta_ij + beta_ji + gamma_ij + gamma_ji - alpha_ij - alpha_ji;
            if (sum < -1e-12) {
                e.flip();
                numFlips++;
                edgeLengthsCM[e] = (l_gamma_ij * l_beta_ji + l_gamma_ji * l_beta_ij) / l_alpha_ij;
            }
        }
    } while (numFlips > 0);
}

FaceData<double> Operators::computeAreas(HalfedgeMesh* mesh, EdgeData<double> &edgeLengths) {
    FaceData<double> faceAreas(mesh);
    for (FacePtr f : mesh->faces()) {
        double l_ij = edgeLengths[f.halfedge().edge()       ];
        double l_jk = edgeLengths[f.halfedge().next().edge()];
        double l_ki = edgeLengths[f.halfedge().prev().edge()];

        double s = (l_ij + l_jk + l_ki) / 2.0;
        faceAreas[f] = sqrt( s * (s - l_ij) * (s - l_jk) * (s - l_ki) );
    }
    return faceAreas;
}


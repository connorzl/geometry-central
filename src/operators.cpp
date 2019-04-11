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


#include "geometrycentral/uniformization.h"

Uniformization::Uniformization(HalfedgeMesh* m, Geometry<Euclidean> *g, VertexData<int> s) : mesh(m), geom(g), singularities(s) {
    cmEdgeLengths = EdgeData<double>(mesh);
    cmAngles = HalfedgeData<double>(mesh);
    cmAreas = FaceData<double>(mesh);
    cmCurvatures = VertexData<double>(mesh);
    cmEdgeVector = HalfedgeData<std::complex<double>>(mesh);
}

void Uniformization::uniformize() {
    if (mesh->eulerCharacteristic() != 2) {
        throw std::logic_error("non genus-0 shape");
    }
    if (mesh->nBoundaryLoops() > 0) {
        //throw std::logic_error ("not implemented for boundary meshes yet");
        uniformizeBoundary();
        computeEdgeVectors();
        return;
    }

    std::cout << "Boundary-less Uniformization..." << std::endl;
    VertexData<size_t> vertexIndices = mesh->getVertexIndices();

    // Ax = b, where A = |V| x |V| laplacian, x = |V| x 1 vertex scaling, b = |V| x 1 curvature diff K - K*
    Eigen::MatrixXd KTarg(mesh->nVertices(),1);
    for (VertexPtr v : mesh->vertices()) {
        size_t index = vertexIndices[v];
        if (singularities[v] == 0) {
            KTarg(index,0) = 0;
        } else if (singularities[v] == 1) {
            KTarg(index,0) = M_PI_2;
        } else if (singularities[v] == -1) {
            KTarg(index,0) = -M_PI_2;
        } else {
            throw std::logic_error("non +1/-1 singular index");
        }
    }
    if (std::abs(KTarg.sum() - mesh->eulerCharacteristic() * 2.0 * M_PI) > 1e-8) {
        throw std::logic_error("Target curvatures do not satisfy Gauss Bonnet");
    }

    // make mesh delaunay
    geom->getEdgeLengths(cmEdgeLengths);    
    Operators::hyperbolicEdgeFlips(mesh, cmEdgeLengths);

    Eigen::SparseMatrix<double> A;
    Eigen::MatrixXd x, K;
    double resid = 0;
    do {
        A = Operators::intrinsicLaplaceMatrix(mesh,cmEdgeLengths);
        K = Operators::intrinsicCurvature(mesh,cmEdgeLengths);

        Eigen::MatrixXd rhs = KTarg - K;
        x = solveSquare<double>(A, rhs);
        x = x.array() - x.mean();

        // update edge lengths
        for (EdgePtr e : mesh->edges()) {
            VertexPtr vi = e.halfedge().vertex();
            VertexPtr vj = e.halfedge().twin().vertex();

            double ui = x(vertexIndices[vi],0);
            double uj = x(vertexIndices[vj],0);
            double s = std::exp( (ui + uj) / 2 );
            cmEdgeLengths[e] *= s;
        }

        // perform hyperbolic edge flips to make sure things are okay
        Operators::hyperbolicEdgeFlips(mesh, cmEdgeLengths);

        // check to see if triangle inequality still holds
        for (FacePtr f : mesh->faces()) {
            HalfedgePtr he = f.halfedge();
            double a = cmEdgeLengths[he.edge()];
            double b = cmEdgeLengths[he.next().edge()];
            double c = cmEdgeLengths[he.prev().edge()];
            if (a > b + c || b > a + c || c > a + b) {
                throw std::logic_error("Triangle Inequality Violated during Uniformization!");
            }
        }

        K = Operators::intrinsicCurvature(mesh,cmEdgeLengths);
        resid = (KTarg - K).array().abs().maxCoeff();
        std::cout << "Resid: " << resid << std::endl;  
    } while (resid > 1e-12);

    // store curvatures for visualization
    for (VertexPtr v : mesh->vertices()) {
        cmCurvatures[v] = K(vertexIndices[v],0);        
    }

    // update CM areas and angles
    cmAreas = Operators::computeAreas(mesh, cmEdgeLengths);
    cmAngles = Operators::computeAngles(mesh, cmEdgeLengths);
    computeEdgeVectors();
    std::cout << "Done!" << std::endl;
}
 
// this needs to include hyperbolic edge flips
void Uniformization::uniformizeBoundary() {
    std::cout << "Boundary Uniformization..." << std::endl;
    // re-index the interior vertices starting from 0
    // we want to omit boundary vertices from the energy and pin u to 0 
    size_t iN = 0;
    VertexData<size_t> vertexIndices = mesh->getVertexIndices();
    VertexData<size_t> interiorVertexIndices(mesh);
    Eigen::Array<bool,Eigen::Dynamic,1> interior(mesh->nVertices(),1);

    for (VertexPtr v : mesh->vertices()) {
        if (!v.isBoundary()) {
            interiorVertexIndices[v] = iN++;
        }
        interior(vertexIndices[v],0) = (!v.isBoundary());
    }
    if (iN != mesh->nInteriorVertices()) {
        throw std::logic_error("error indexing interior vertices!");
    }
    
    // Ax = b, where A = |V| x |V| laplacian, x = |V| x 1 vertex scaling, b = |V| x 1 curvature diff K - K*
    Eigen::MatrixXd KTarg(mesh->nInteriorVertices(),1);
    Eigen::MatrixXd K(mesh->nInteriorVertices(),1);
    for (VertexPtr v : mesh->vertices()) {
        if (v.isBoundary()) continue;
        size_t index = interiorVertexIndices[v];
        if (singularities[v] == 0) {
            KTarg(index,0) = 0;
        } else if (singularities[v] == 1) {
            KTarg(index,0) = M_PI_2;
        } else if (singularities[v] == -1) {
            KTarg(index,0) = -M_PI_2;
        } else {
            throw std::logic_error("non +1/-1 singularity");
        }
    }

    // make mesh delaunay
    geom->getEdgeLengths(cmEdgeLengths);
    Operators::hyperbolicEdgeFlips(mesh, cmEdgeLengths);
    
    double resid = 0;
    do {
        Eigen::SparseMatrix<double> A_all = Operators::intrinsicLaplaceMatrix(mesh,cmEdgeLengths);
        BlockDecompositionResult<double> B = blockDecomposeSquare(A_all, interior);

        Eigen::MatrixXd K_all = Operators::intrinsicCurvature(mesh,cmEdgeLengths);
        for (VertexPtr v : mesh->vertices()) {
            if (v.isBoundary()) continue;
            K(interiorVertexIndices[v],0) = K_all(vertexIndices[v],0);
        }

        Eigen::MatrixXd rhs = KTarg - K;
        Eigen::MatrixXd x = solveSquare<double>(B.AA, rhs);

        // update edge lengths
        for (EdgePtr e : mesh->edges()) {
            VertexPtr vi = e.halfedge().vertex();
            VertexPtr vj = e.halfedge().twin().vertex();

            double ui, uj;
            if (vi.isBoundary()) {
                ui = 0;
            } else {
                ui = x(interiorVertexIndices[vi],0);
            }
            if (vj.isBoundary()) {
                uj = 0;
            } else {
                uj = x(interiorVertexIndices[vj],0);
            }
            double s = std::exp( (ui + uj) / 2 );
            cmEdgeLengths[e] *= s;
        }
        
        Operators::hyperbolicEdgeFlips(mesh, cmEdgeLengths);

        // check to see if triangle inequality still holds
        for (FacePtr f : mesh->faces()) {
            HalfedgePtr he = f.halfedge();
            double a = cmEdgeLengths[he.edge()];
            double b = cmEdgeLengths[he.next().edge()];
            double c = cmEdgeLengths[he.prev().edge()];
            if (a > b + c || b > a + c || c > a + b) {
                throw std::logic_error("Triangle Inequality Violated during Uniformization!");
            }
        }

        K_all = Operators::intrinsicCurvature(mesh,cmEdgeLengths);
        for (VertexPtr v : mesh->vertices()) {
            if (v.isBoundary()) continue;
            K(interiorVertexIndices[v],0) = K_all(vertexIndices[v],0);
        }
        resid = (KTarg - K).array().abs().maxCoeff();
        std::cout << "Resid: " << resid << std::endl;  
    } while (resid > 1e-12);

    K = Operators::intrinsicCurvature(mesh, cmEdgeLengths);
    // store curvatures for visualization
    for (VertexPtr v : mesh->vertices()) {
        cmCurvatures[v] = K(vertexIndices[v],0);        
    }

    // update CM areas and angles
    cmAreas = Operators::computeAreas(mesh, cmEdgeLengths);
    cmAngles = Operators::computeAngles(mesh, cmEdgeLengths);
    std::cout << "Done!" << std::endl;
}

void Uniformization::computeEdgeVectors() {
    for (FacePtr f : mesh->faces()) {
        // Gather elements
        HalfedgePtr he = f.halfedge();
        double l_jl = cmEdgeLengths[he.edge()];
        double l_ij = cmEdgeLengths[he.prev().edge()];
        double theta_ijl = cmAngles[he.next()]; // radians

        // Place first vertex at (0,0)
        Vector2 j = Vector2{0,0};
        Vector2 l = Vector2{l_jl,0};
        Vector2 i = Vector2{cos(theta_ijl) * l_ij, sin(theta_ijl) * l_ij}; 

        Vector2 jl = l - j;
        Vector2 li = i - l;
        Vector2 ij = j - i;

        cmEdgeVector[he] = std::complex<double>(jl.x,jl.y);
        cmEdgeVector[he.next()] = std::complex<double>(li.x,li.y);
        cmEdgeVector[he.prev()] = std::complex<double>(ij.x,ij.y);
    }
}
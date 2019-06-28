#pragma once

#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/geometry.h"
#include "geometrycentral/operators.h"
#include "geometrycentral/linear_solvers.h"

using namespace geometrycentral;

class Uniformization {
public:
    Uniformization() {};
    Uniformization(HalfedgeMesh* m, Geometry<Euclidean> *g, VertexData<int> s);
    void uniformize();

    EdgeData<double> cmEdgeLengths;
    HalfedgeData<double> cmAngles;
    FaceData<double> cmAreas;
    VertexData<double> cmCurvatures;

    HalfedgeData<std::complex<double>> cmEdgeVector;

private:
    HalfedgeMesh* mesh;
    Geometry<Euclidean>* geom;
    VertexData<int> singularities;
    
    void uniformizeBoundary();

    void computeEdgeVectors();
};
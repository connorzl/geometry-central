#pragma once

#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/geometry.h"

using namespace geometrycentral;

class Operators {
public:
    static Eigen::SparseMatrix<double> laplaceMatrix(HalfedgeMesh* m, Geometry<Euclidean>* g);
    static Eigen::SparseMatrix<double> intrinsicLaplaceMatrix(HalfedgeMesh* m, EdgeData<double> &edgeLengths);
    static Eigen::MatrixXd intrinsicCurvature(HalfedgeMesh* m, EdgeData<double> &edgeLengths);
    static Eigen::SparseMatrix<double> buildExteriorDerivative0Form(HalfedgeMesh* m, Geometry<Euclidean>* g);
    static Eigen::SparseMatrix<double> buildHodgeStar1Form(HalfedgeMesh* m, Geometry<Euclidean>* g);
    static HalfedgeData<double> computeAngles(HalfedgeMesh *m, EdgeData<double> &edgeLengths);
    static FaceData<double> computeAreas(HalfedgeMesh *m, EdgeData<double> &edgeLengths);
    static void hyperbolicEdgeFlips(HalfedgeMesh* mesh, EdgeData<double> &edgeLengthsCM);
};
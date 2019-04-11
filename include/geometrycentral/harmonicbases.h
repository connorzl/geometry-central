#pragma once

#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/geometry.h"
#include "geometrycentral/polygon_soup_mesh.h"
#include "geometrycentral/operators.h"

using namespace geometrycentral;

class HarmonicBases {
    public:
        HarmonicBases(HalfedgeMesh* m, Geometry<Euclidean>* g);
        std::vector<Eigen::MatrixXd> compute();
        FaceData<std::complex<double>> visualize();

        std::vector<std::vector<HalfedgePtr>> generators;
        std::vector<Eigen::MatrixXd> bases;

    private:
        HalfedgeMesh* mesh;
        Geometry<Euclidean>* geom;

        // quantities for tree-cotree
        std::map<VertexPtr,VertexPtr> vertexParent;
        std::map<FacePtr,FacePtr> faceParent;
        

        // helper functions for tree-cotree
        void treeCotree();
        void buildPrimalSpanningTree();
        bool inPrimalSpanningTree(HalfedgePtr he);
        void buildDualSpanningCotree();
        bool inDualSpanningCotree(HalfedgePtr he);

        // quantities for harmonic bases
        HalfedgeData<int> orientation;
        Eigen::SparseMatrix<double> A;
        Eigen::SparseMatrix<double> d0;
        Eigen::SparseMatrix<double> d0T;
        Eigen::SparseMatrix<double> hodge1;

        // helper functions for harmonic bases
        Eigen::MatrixXd buildClosedPrimalOneForm(std::vector<HalfedgePtr> generator);
        Eigen::MatrixXd computeExactComponent(Eigen::MatrixXd omega);
};
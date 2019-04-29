#pragma once

#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/geometry.h"
#include "geometrycentral/polygon_soup_mesh.h"

using namespace geometrycentral;

class Stripes {
    public:
        Stripes(HalfedgeMesh* m, Geometry<Euclidean>* g);
        VertexData<std::complex<double>> computeField();
        FaceData<int> computeSingularities();
        void edgeData();
        void computeStripes();

    private:
        HalfedgeMesh* mesh;
        Geometry<Euclidean>* geom;

        int nPowerIterations = 20;
        double eps = std::pow(10.,-8.);
        double n = 2;
        double lambda = 130;

        // quantities
        HalfedgeData<double> phi; // this is the theta referred to in the stripes paper
        HalfedgeData<std::complex<double>> r;
        VertexData<std::complex<double>> field;
        FaceData<int> singularities;
        EdgeData<int> branchCover;
        EdgeData<double> omega;
        
        // helpers
        void setup();
        Eigen::SparseMatrix<std::complex<double>> assembleM();
        Eigen::SparseMatrix<std::complex<double>> assembleA();
        Eigen::MatrixXcd principalEigenvector(Eigen::SparseMatrix<std::complex<double>> A, Eigen::SparseMatrix<std::complex<double>> B);
        Eigen::SparseMatrix<double> EnergyMatrix();
        Eigen::SparseMatrix<double> MassMatrix();
        Eigen::MatrixXcd principalEigenvector(Eigen::SparseMatrix<double> A, Eigen::SparseMatrix<double> B);

        // visualization
        //Eigen::SparseMatrix<std::complex<double>> buildLaplacian();
};
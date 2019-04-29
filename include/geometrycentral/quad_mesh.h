#pragma once

#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/geometry.h"
#include "geometrycentral/polygon_soup_mesh.h"
#include "geometrycentral/operators.h"
#include "geometrycentral/BHalfedgemesh.h"

using namespace geometrycentral;

class QuadMesh {
    public:
        QuadMesh(HalfedgeMesh* m, Geometry<Euclidean>* g);
        FaceData<std::complex<double>> computeCrossField();
        VertexData<int> computeSingularities();
        HalfedgeData<int> computeBranchCover(bool improve=true);
        
        // cone metric
        VertexData<double> uniformize();
        FaceData<std::complex<double>> computeCrossFieldCM();
        std::vector<FaceData<std::complex<double>>> computeCrossFieldCMBranchCover();
        
        // stripes
        std::vector<VertexData<double>> computeStripes();
        std::vector<EdgeData<double>> omega;
        std::vector<VertexData<double>>  visualize();

    private:
        HalfedgeMesh* mesh;
        Geometry<Euclidean>* geom;

        int nPowerIterations = 20;
        double eps = std::pow(10.,-8.);
        int n = 4;

        // smoothest field quantities
        HalfedgeData<std::complex<double>> theta;
        HalfedgeData<std::complex<double>> r;
        FaceData<std::complex<double>> field;

        // smoothest field helpers
        void setup();
        Eigen::SparseMatrix<std::complex<double>> assembleM();
        Eigen::SparseMatrix<std::complex<double>> assembleA();
        void computeSmoothestField(Eigen::SparseMatrix<std::complex<double>> M, Eigen::SparseMatrix<std::complex<double>> A);

        VertexData<int> singularities;
        HalfedgeData<int> eta;
        int numSingularities = 0;

        // uniformization results
        EdgeData<double> edgeLengthsCM;
        HalfedgeData<std::complex<double>> thetaCM;
        HalfedgeData<std::complex<double>> rCM;
        HalfedgeData<double> cmAngles;
        
        // cross field
        FaceData<std::complex<double>> fieldCM;    // field computed on original surface
        std::vector<FaceData<std::complex<double>>> branchCoverFields; // field computed on branch cover
        BranchCoverTopology BC;
        double scale = 1;

        // stripes
        void computeOmega();
        Eigen::SparseMatrix<double> energyMatrix();
        Eigen::SparseMatrix<double> massMatrix();
        std::vector<VertexData<size_t>> BVertexIndices;

        // visualization
        std::vector<VertexData<std::complex<double>>> psi;
        std::vector<VertexData<double>> coords;

        // helpers
        void setupCM();
};
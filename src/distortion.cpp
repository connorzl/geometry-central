#include "geometrycentral/distortion.h"

using namespace geometrycentral;
using std::cout;
using std::endl;

std::vector<double> Distortion::distortion;

double computeAreaScaling(const std::vector<Vector3>& p, const std::vector<Vector3>& q) {
    // Compute edge vectors and areas
    Vector3 u1 = p[1] - p[0]; 
    Vector3 u2 = p[2] - p[0]; 
    double Area = norm(cross(u1, u2));

    Vector3 v1 = q[1] - q[0];
    Vector3 v2 = q[2] - q[0];
    double area = norm(cross(v1, v2));

    return log(area / Area);
}

Vector3 Distortion::computeAreaScaling(HalfedgeMesh* mesh, Geometry<Euclidean>* geom) {
    double totalU = 0.0;
    double minU = std::numeric_limits<double>::infinity();
    double maxU = -std::numeric_limits<double>::infinity();
    double totalArea = 0.0;
    distortion.resize(mesh->nFaces());
    std::vector<Vector3> p(3), q(3);

    for (size_t i = 0; i < mesh->nFaces(); i++) {
        FacePtr f = mesh->face(i);
        HalfedgePtr he = f.halfedge();
        int j = 0;
        do {
            p[j] = geom->position(he.vertex());
            Vector2 uv = geom->paramCoords[he.next()];
            q[j] = Vector3{uv.x, uv.y, 0.0};
            j++;

            he = he.next();
        } while (he != f.halfedge());

        double u = ::computeAreaScaling(p, q);
        double area = geom->area(f);
        totalU += u * area;
        totalArea += area;

        // Clamp to range [0.0, 1.0]
        distortion[i] = std::max(0.0, std::min(1.0, u));
        minU = std::min(minU, u);
        maxU = std::max(maxU, u);
    }

    return Vector3 { minU, maxU, totalU / totalArea };
}

double computeQuasiConformalerror(std::vector<Vector3>& p, std::vector<Vector3>& q) {
    // Compute edge vectors
    Vector3 u1 = p[1] - p[0];
    Vector3 u2 = p[2] - p[0];

    Vector3 v1 = q[1] - q[0];
    Vector3 v2 = q[2] - q[0];

    // Compute orthonormal bases
    Vector3 e1 = u1;
    e1.normalize();
    Vector3 e2 = u2 - dot(u2, e1) * e1;
    e2.normalize();

    Vector3 f1 = v1;
    f1.normalize();
    Vector3 f2 = v2 - dot(v2, f1) * f1;
    f2.normalize();

    // Project onto bases
    p[0] = Vector3{0, 0, 0};
    p[1] = Vector3{dot(u1, e1), dot(u1, e2), 0};
    p[2] = Vector3{dot(u2, e1), dot(u2, e2), 0};

    q[0] = Vector3{0, 0, 0};
    q[1] = Vector3{dot(v1, f1), dot(v1, f2), 0};
    q[2] = Vector3{dot(v2, f1), dot(v2, f2), 0};

    double A = 2.0 * norm(cross(u1, u2));

    Vector3 Ss = (q[0]*(p[1].y - p[2].y) + q[1]*(p[2].y - p[0].y) + q[2]*(p[0].y - p[1].y)) / A;
    Vector3 St = (q[0]*(p[2].x - p[1].x) + q[1]*(p[0].x - p[2].x) + q[2]*(p[1].x - p[0].x)) / A;
    double a = dot(Ss, Ss);
    double b = dot(Ss, St);
    double c = dot(St, St);
    
    double det = sqrt(pow(a-c, 2) + 4.0*b*b);
    double Gamma = sqrt(0.5*(a + c + det));
    double gamma = sqrt(0.5*(a + c - det));

    if (Gamma < gamma) std::swap(Gamma, gamma);
    return Gamma/gamma;
}

Vector3 Distortion::computeQuasiConformalError(HalfedgeMesh* mesh, Geometry<Euclidean>* geom) {
    double totalQc = 0.0;
    double minQc = std::numeric_limits<double>::infinity();
    double maxQc = -std::numeric_limits<double>::infinity();
    double totalArea = 0.0;
    distortion.resize(mesh->nFaces());
    std::vector<Vector3> p(3), q(3);

    for (size_t i = 0; i < mesh->nFaces(); i++) {
        FacePtr f = mesh->face(i);
        HalfedgePtr he = f.halfedge();
        int j = 0;
        do {
            p[j] = geom->position(he.vertex());
            Vector2 uv = geom->paramCoords[he.next()];
            q[j] = Vector3 {uv.x, uv.y, 0.0};
            j++;

            he = he.next();
        } while (he != f.halfedge());

        double qc = ::computeQuasiConformalerror(p, q);
        double area = geom->area(f);
        totalQc += qc * area;
        totalArea += area;

        // Clamp distortion to [1, 1.5]
        distortion[i] = std::max(1.0, std::min(1.5, qc));
        maxQc = std::max(maxQc, qc);
        minQc = std::min(minQc, qc);
    }

    return Vector3 { minQc, maxQc, totalQc / totalArea};
}
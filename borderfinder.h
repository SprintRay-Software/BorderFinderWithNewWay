#ifndef BORDERFINDER_H
#define BORDERFINDER_H

#include "borderfinder_global.h"
#include "Clipper.h"
#include <qpoint.h>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>

using namespace std;

typedef OpenMesh::PolyMesh_ArrayKernelT<> Mesh;
typedef Mesh::Point MyPoint;
typedef std::vector<MyPoint> Boundary;
//Mesh myMesh1;

constexpr int precision = 3;

class BORDERFINDERSHARED_EXPORT BorderFinder
{

public:
    BorderFinder();
    BorderFinder(Mesh myMesh);
    bool startFinder(QString src_path, double renderScale, bool isRead, QString modelFilePath);
    bool cloneModelAndProcess(string inputPath,double scale, QString modelFilePath);
    bool calcOff(int &offX, int &offY,double scale);
    bool write_mesh(const QString& filename);
    void ramerDouglasPeucker(const vector<QPointF> &pointList, double epsilon, vector<QPointF> &out);
    double perpendicularDistance(const QPointF &pt, const QPointF lineStart, const QPointF lineEnd);
    bool drawBorderPointsToImage(QVector<QPointF> pointList,QString imagePath);
    Mesh ax(string inputPath);

public:
    int offX;
    int offY;
    Mesh myMesh;
    Mesh myMeshBuff;
    bool isReload = false;
};

#endif // BORDERFINDER_H

#include <QCoreApplication>
#include <qdir.h>
#include "borderfinder.h"
#include <QTextStream>
#include <QDebug>
#include <QPainter>
#include <QTextCodec>
#include "borderfinderConstants.h"
#include <QTextCodec>
#include <QtMath>
#include "grahamscan.h"


BorderFinder::BorderFinder()
{
}

BorderFinder::BorderFinder(Mesh myMesh)
{
    this->myMesh = myMesh;
}

bool BorderFinder::startFinder(QString src_path, double renderScale, bool isRead, QString modelFilePath)
{
    QTextCodec *code = QTextCodec::codecForName("GBK");
    myMesh = myMeshBuff;
    if(myMeshBuff.vertices_empty())
    {
        qDebug()<<"myMesh is empty!";
    }
    if(isRead)
    {
        myMesh.clear();
        myMesh.request_face_normals();
        string str = code->fromUnicode(src_path).data();
        if (!OpenMesh::IO::read_mesh(myMesh, str))
        {
            std::cout << "Fail to read mesh!" << std::endl;
            return false;
        }
        myMeshBuff = myMesh;
    }
    //QString tempPath1 = QCoreApplication::applicationDirPath()+"/Borderfinder/C1.stl";
    string tempPath1 = code->fromUnicode(QCoreApplication::applicationDirPath()+"/Borderfinder/C1.stl").data();
    OpenMesh::IO::write_mesh(myMesh, tempPath1, OpenMesh::IO::Options::Binary);
    OpenMesh::IO::read_mesh(myMesh, tempPath1);

    //myMesh.halfedges();
    //myMesh.set_halfedge_handle();
    myMesh.request_vertex_status();
    myMesh.request_halfedge_status();
    myMesh.request_edge_status();
    myMesh.request_face_status();
    myMesh.request_face_normals();
    myMesh.update_normals();
    const MyPoint standard_normal(0, 0, 1);
    for (auto f_it = myMesh.faces_begin(); f_it != myMesh.faces_end(); f_it++)
    {
//            for(auto fv_it = myMesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
//            {
//                myMesh.point (*fv_it).data()[0] = abs(myMesh.point (*fv_it).data()[0]);
//                myMesh.point (*fv_it).data()[1] = abs(myMesh.point (*fv_it).data()[1]);
//                myMesh.point (*fv_it).data()[2] = abs(myMesh.point (*fv_it).data()[2]);
//            }
        auto p =myMesh.normal(*f_it);
        double x = p[0];
        double y = p[1];
        double z = p[2];
        if ((x*standard_normal[0] + y * standard_normal[1] + z * standard_normal[2]) < 0)
        {
            myMesh.delete_face(*f_it);
        }
    }
    myMesh.garbage_collection();
    std::vector<Boundary> boundaries;

    int j = 0;
    for (auto he_it = myMesh.halfedges_begin(); he_it != myMesh.halfedges_end(); he_it++)
    {
        j++;
        if (myMesh.is_boundary(*he_it) && !myMesh.status(*he_it).locked())
        {
            Boundary boundary;
            Mesh::HalfedgeHandle heh_init = *he_it;
            Mesh::HalfedgeHandle heh_curr = *he_it;
            do
            {
                myMesh.status(heh_curr).set_locked(true);
                boundary.push_back(myMesh.point(myMesh.from_vertex_handle(heh_curr)));
                heh_curr = myMesh.next_halfedge_handle(heh_curr);
            } while (heh_curr != heh_init);
            boundaries.push_back(std::move(boundary));
        }
    }
    ClipperLib::Clipper clipper;
    ClipperLib::Paths src_paths;
    ClipperLib::Paths dst_paths;
    ClipperLib::ClipperOffset  clipperoff;
    const double scale = std::pow(10, precision);
    for (auto& boundary : boundaries)
    {
        ClipperLib::Path dst_path;
        for (const auto& point : boundary)
        {
            dst_path.push_back(ClipperLib::IntPoint(static_cast<ClipperLib::cInt>(point[0] * scale), static_cast<ClipperLib::cInt>(point[1] * scale)));
        }
        if (ClipperLib::Orientation(dst_path))
            continue;
        src_paths.push_back(std::move(dst_path));
    }
    //ClipperLib::SimplifyPolygons(src_paths,src_paths,ClipperLib::pftNonZero);
    clipper.AddPaths(src_paths, ClipperLib::PolyType::ptSubject, true);
    clipper.Execute(ClipperLib::ClipType::ctUnion, dst_paths, ClipperLib::PolyFillType::pftNonZero);
    ClipperLib::CleanPolygons(dst_paths,10);
    clipperoff.AddPaths(dst_paths,ClipperLib::jtSquare,ClipperLib::etClosedPolygon);
    clipperoff.Execute(dst_paths,0 * scale);
    ClipperLib::SimplifyPolygons(dst_paths,dst_paths, ClipperLib::PolyFillType::pftNonZero);
    boundaries.clear();
    for (const auto& dst_path : dst_paths)
    {
        if (dst_path.size() == 0)
            continue;
        Boundary boundary;
        for (const auto& point : dst_path)
        {
            boundary.push_back(MyPoint(point.X / scale, point.Y / scale, 0.0));
        }
        //boundary.push_back(boundary[0]);
        boundaries.push_back(std::move(boundary));
    }

//    ClipperLib::Paths src_paths1;
//    ClipperLib::Paths dst_paths1;
//    for (auto& boundary : boundaries)
//    {
//        ClipperLib::Path dst_path;
//        for (const auto& point : boundary)
//        {
//            dst_path.push_back(ClipperLib::IntPoint(static_cast<ClipperLib::cInt>(point[0] * scale), static_cast<ClipperLib::cInt>(point[1] * scale)));
//        }
//        src_paths1.push_back(std::move(dst_path));
//    }
//    //ClipperLib::SimplifyPolygons(src_paths,src_paths,ClipperLib::pftNonZero);
//    clipper.AddPaths(src_paths1, ClipperLib::PolyType::ptSubject, true);
//    clipper.Execute(ClipperLib::ClipType::ctUnion, dst_paths1, ClipperLib::PolyFillType::pftEvenOdd);

//    boundaries.clear();
//    for (const auto& dst_path : dst_paths1)
//    {
//        if (!ClipperLib::Orientation(dst_path))
//            continue;
//        Boundary boundary;
//        for (const auto& point : dst_path)
//        {
//            boundary.push_back(MyPoint(point.X / scale, point.Y / scale, 0.0));
//        }

//        boundaries.push_back(std::move(boundary));
//    }

    QFileInfo fileInfo(src_path);
    QString fileName = fileInfo.baseName();
    //QString fileName = src_path.right(src_path.size()-src_path.lastIndexOf('.')-1);
    if(fileName != modelFilePath)
    {
        fileName = modelFilePath;
    }
    qDebug()<<"startFinder->writepath: "<<QCoreApplication::applicationDirPath() + + "/Borderfinder/" + fileName + ".txt", boundaries;

    QString filepng(QCoreApplication::applicationDirPath() + + "/Borderfinder/" + fileName + ".png");
    QVector<QPointF> pointList1;
    for (const auto& boundary : boundaries)
    {
        for (const auto& point : boundary)
        {           
            QPointF tempPoint(point[0]*(1.0/renderScale),point[1]*(1.0/renderScale));
            pointList1.push_back(tempPoint);
        }
        //break;
    }
    for (auto point : pointList1)
    {
        //qDebug()<<"drawBorderPointsToImage.pointListOut: "<<point;
    }

    if(!drawBorderPointsToImage(pointList1,filepng))
    {
        return false;
    }

    QFile file(QCoreApplication::applicationDirPath() + + "/Borderfinder/" + fileName + ".txt");
    if (file.open(QFile::WriteOnly)) {
        QTextStream out(&file);

        vector<QPointF> clonePointList;
        vector<QPointF> pointList;
        for (const auto& boundary : boundaries)
        {
            for (const auto& point : boundary)
            {
                QPointF tempPoint(point[0]*(1.0/renderScale),point[1]*(1.0/renderScale));
                clonePointList.push_back(tempPoint);
            }
            //break;
        }

        if (clonePointList.size() <= 0)
        {
            return false;
        }

        ramerDouglasPeucker(clonePointList, 15.0, pointList);

        double minX = pointList[0].x(), maxX = pointList[0].x(), minY = pointList[0].y(), maxY = pointList[0].y();
        for (QPointF point : pointList)
        {
            if (point.x() <= minX)
            {
                minX = point.x();
            }
            if (point.y() <= minY)
            {
                minY = point.y();
            }
            if (point.x() >= maxX)
            {
                maxX = point.x();
            }
            if (point.y() >= maxY)
            {
                maxY = point.y();
            }
        }

        double coorW = (maxX + minX) / 2.0;
        double coorH = (maxY + minY) / 2.0;

        double y_offset = 1080 - 2 * coorH;

        out << pointList.size();
        QVector<QPointF> pointList3;
        for (QPointF point : pointList) {
             out << ' ' << '[' << point.x() << ',' << 2*coorH - point.y() << ']';
             pointList3.push_back(QPointF(point.x(), 2*coorH - point.y()));
        }
        out << endl;
        file.close();

//        if(!drawBorderPointsToImage(pointList3,filepng))
//        {
//            return false;
//        }
    }
//    myMesh.clear();
//    myMeshBuff.clear();
    return true;
}

Mesh BorderFinder::ax(string inputPath)
{
    Mesh mesh1;
    QTextCodec *code = QTextCodec::codecForName("GBK");
    string str = code->fromUnicode(QString::fromStdString(inputPath)).data();
    if (!OpenMesh::IO::read_mesh(myMesh, str))
    {
        std::cout << "Failed to read mesh!" << std::endl;
    }
    return mesh1;
}

bool BorderFinder::cloneModelAndProcess(string inputPath,double scale, QString modelFilePath)
{
    QTextCodec *code = QTextCodec::codecForName("GBK");
    string str = code->fromUnicode(QString::fromStdString(inputPath)).data();
    if (!OpenMesh::IO::read_mesh(myMesh, str))
    {
        std::cout << str << std::endl;
        std::cout << "Failed to read mesh!" << std::endl;
        return false;
    }

    offX = 0, offY = 0;
    bool isNormal = calcOff(offX, offY, scale);
    if(!isNormal)
    {
        return isNormal;
    }

    QString generatePath1 = QCoreApplication::applicationDirPath() + BORDERFINDERCONSTANTS::CLONE_STL_PATH + "/d.stl";
    auto localPath1 = generatePath1.toLocal8Bit();
    str = code->fromUnicode(localPath1).data();
    OpenMesh::IO::write_mesh(myMesh, str, OpenMesh::IO::Options::Binary);

    if(offX != 0 || offY != 0)
    {
        QTextCodec *code = QTextCodec::codecForName("GBK");
        string str = code->fromUnicode(QString::fromStdString(inputPath)).data();
        QString qinputPath = QString(QString::fromLocal8Bit(str.c_str()));
        QDir dir(QCoreApplication::applicationDirPath() + BORDERFINDERCONSTANTS::SQUASHED_IMG_PATH);
        QFileInfo fileInfo(qinputPath);
        QString fileName = fileInfo.baseName();
        //QString fileName = qinputPath.right(qinputPath.size()-qinputPath.lastIndexOf('.')-1);
        if(fileName != modelFilePath)
        {
            fileName = modelFilePath;
        }
        bool exist = dir.exists(fileName);
        if (!exist)
        {
            dir.mkdir(fileName);
        }
        QString generatePath = QCoreApplication::applicationDirPath() + BORDERFINDERCONSTANTS::CLONE_STL_PATH + "/" + fileName + "/" + fileName + ".stl";
        //auto localPath = generatePath.toLocal8Bit();
        //qDebug()<<"read_mesh11111: "<<QString::fromStdString(inputPath);
        //qDebug()<<"write_mesh11111: "<<localPath;
        str = code->fromUnicode(generatePath).data();
        if (!OpenMesh::IO::write_mesh(myMesh, str, OpenMesh::IO::Options::Binary))
        {
            std::cerr << "Cannot write mesh to file ' .stl ' " << std::endl;
            return false;
        }
        startFinder(generatePath, scale, true, modelFilePath);
        QString path = QCoreApplication::applicationDirPath() + BORDERFINDERCONSTANTS::CLONE_STL_PATH + "/" + fileName;
        dir = QDir(path.toLocal8Bit());
        foreach(QFileInfo mfi, dir.entryInfoList())
        {
            if (mfi.isFile() && (mfi.suffix() == "stl1"))
            {
                dir.remove(mfi.fileName().toLocal8Bit());
            }
        }
        isReload = true;
        return false;
    }
    else
    {
        startFinder(QString::fromStdString(inputPath), scale, true, modelFilePath);
    }
    return isNormal;
}
bool BorderFinder::calcOff(int &offX, int &offY, double scale)
{
    if(myMesh.vertices_empty())
    {
        return false;
    }
    float minX = myMesh.point(myMesh.vertices_begin()).data()[0], maxX = myMesh.point(myMesh.vertices_begin()).data()[0],
        minY = myMesh.point(myMesh.vertices_begin()).data()[1], maxY = myMesh.point(myMesh.vertices_begin()).data()[1];
    for (auto vec_it = myMesh.vertices_begin(); vec_it != myMesh.vertices_end(); vec_it++)
    {
        if (myMesh.point(*vec_it).data()[0] <= minX)
        {
            minX = myMesh.point(*vec_it).data()[0];
        }
        if (myMesh.point(*vec_it).data()[0] >= maxX)
        {
            maxX = myMesh.point(*vec_it).data()[0];
        }
        if (myMesh.point(*vec_it).data()[1] <= minY)
        {
            minY = myMesh.point(*vec_it).data()[1];
        }
        if (myMesh.point(*vec_it).data()[1] >= maxY)
        {
            maxY = myMesh.point(*vec_it).data()[1];
        }
    }
    double length = 1920 * scale;
    double width = 1080 * scale;
    if (minX < 0)
    {
        offX = -((int)minX - 1);
    }
    if (maxX > length)
    {
        offX = -(int)(maxX - length + 1);
    }
    if (minY < 0)
    {
        offY = -((int)minY - 1);
    }
    if (maxY > width)
    {
        offY = -(int)(maxY - width + 1);
    }
    std::cout << "offX = " <<offX<< std::endl;
    std::cout << "offY = " <<offY<< std::endl;
    if (offX != 0 || offY != 0)
    {
        for (auto vec_it = myMesh.vertices_begin(); vec_it != myMesh.vertices_end(); vec_it++)
        {
            myMesh.point(vec_it).data()[0] = myMesh.point(vec_it).data()[0] + offX;
            myMesh.point(vec_it).data()[1] = myMesh.point(vec_it).data()[1] + offY;
        }
    }
    return true;
}


void BorderFinder::ramerDouglasPeucker(const vector<QPointF> &pointList, double epsilon, vector<QPointF> &out)
{
    if (pointList.size() < 2)
    {
        return;
    }

    // Find the point with the maximum distance from line between start and end
    double dmax = 0.0;
    size_t index = 0;
    size_t end = pointList.size() - 1;
    for (size_t i = 1; i < end; i++)
    {
        double d = perpendicularDistance(pointList[i], pointList[0], pointList[end]);
        if (d > dmax)
        {
            index = i;
            dmax = d;
        }
    }

    // If max distance is greater than epsilon, recursively simplify
    if (dmax > epsilon)
    {
        // Recursive call
        vector<QPointF> recResults1;
        vector<QPointF> recResults2;
        vector<QPointF> firstLine(pointList.begin(), pointList.begin() + index + 1);
        vector<QPointF> lastLine(pointList.begin() + index, pointList.end());
        ramerDouglasPeucker(firstLine, epsilon, recResults1);
        ramerDouglasPeucker(lastLine, epsilon, recResults2);

        // Build the result list
        out.assign(recResults1.begin(), recResults1.end() - 1);
        out.insert(out.end(), recResults2.begin(), recResults2.end());
        if (out.size() < 2)
        {
            return;
        }
    }
    else
    {
        //Just return start and end points
        out.clear();
        out.push_back(pointList[0]);
        out.push_back(pointList[end]);
    }
}

double BorderFinder::perpendicularDistance(const QPointF &pt, const QPointF lineStart, const QPointF lineEnd)
{
    double dx = lineEnd.x() - lineStart.x();
    double dy = lineEnd.y() - lineStart.y();

    //Normalise
    double mag = pow(pow(dx, 2.0) + pow(dy, 2.0), 0.5);
    if (mag > 0.0)
    {
        dx /= mag; dy /= mag;
    }

    double pvx = pt.x() - lineStart.x();
    double pvy = pt.y() - lineStart.y();

    //Get dot product (project pv onto normalized direction)
    double pvdot = dx * pvx + dy * pvy;

    //Scale line direction vector
    double dsx = pvdot * dx;
    double dsy = pvdot * dy;

    //Subtract this from pv
    double ax = pvx - dsx;
    double ay = pvy - dsy;

    return pow(pow(ax, 2.0) + pow(ay, 2.0), 0.5);
}

bool BorderFinder::drawBorderPointsToImage(QVector<QPointF> pointList,QString imagePath, bool transform)
{
    QImage borderedImg(1920,1080, QImage::Format_RGB32);
    //borderedImg.
    vector<QPointF> clonePointList;
    vector<QPointF> pointListOut;
    for (QPointF point : pointList)
    {
        clonePointList.push_back(point);
    }

    double minX = pointList[0].x(), maxX = pointList[0].x(), minY = pointList[0].y(), maxY = pointList[0].y();
    for (QPointF point : pointList)
    {
        if (point.x() <= minX)
        {
            minX = point.x();
        }
        if (point.y() <= minY)
        {
            minY = point.y();
        }
        if (point.x() >= maxX)
        {
            maxX = point.x();
        }
        if (point.y() >= maxY)
        {
            maxY = point.y();
        }
    }

    double coorW = (maxX + minX) / 2.0;
    double coorH = (maxY + minY) / 2.0;

    double y_offset = 1080 - 2 * coorH;

    ramerDouglasPeucker(clonePointList, 15.0, pointListOut);

    QPainter p(&borderedImg);
    p.setPen(QColor(255, 255, 255));
    int length = pointListOut.size();
    for (int i = 0; i < length; i++)
    {
        if(!transform)
            p.drawLine(pointListOut[i%length].x(), 2*coorH - pointListOut[i%length].y(), pointListOut[(i + 1) % length].x(), 2*coorH - pointListOut[(i + 1) % length].y());
        else
            p.drawLine(pointListOut[i%length].x(), pointListOut[i%length].y(), pointListOut[(i + 1) % length].x(), pointListOut[(i + 1) % length].y());
    }

    // will return false if error(point.y() + y_offset)
    return borderedImg.save(imagePath);
}

int BorderFinder::ComputeByGrahamScan(QVector<QPointF> pOriginalPos, const int nOriginalCount, int iConvexIndices[], int &nConvexPointCount)
{
    const int MAX_POINT_NUMBER = 20;
    if (nOriginalCount < 1 || nOriginalCount > MAX_POINT_NUMBER)
    {
        return -1;
    }

    int iStartPointIndex = 0;

    for (int i=1; i<nOriginalCount; i++)
    {
        if (pOriginalPos[i].y() < pOriginalPos[iStartPointIndex].y())
        {
            iStartPointIndex = i;
        }
        else if (pOriginalPos[i].y() == pOriginalPos[iStartPointIndex].y())
        {
            if (pOriginalPos[i].x() < pOriginalPos[iStartPointIndex].x())
            {
                iStartPointIndex = i;
            }
        }
    }

    float fCrossResult[MAX_POINT_NUMBER] = {0};
    int iCrossIndices[MAX_POINT_NUMBER] = {0};

    for (int i=0; i<nOriginalCount; i++)
    {
        int iReferVector[2] = {1, 0};
        int iCurrentVector[2] = { pOriginalPos[i].x() - pOriginalPos[iStartPointIndex].x(), pOriginalPos[i].y() - pOriginalPos[iStartPointIndex].y()};

        fCrossResult[i] = qAtan((float)iCurrentVector[1]/(float)iCurrentVector[0]);
        //fCrossResult[i] = Actan_0_2PI((float)iCurrentVector[1], (float)iCurrentVector[0]);

        iCrossIndices[i] = i;
    }

    //   说明：极角排序 [10/21/2016 ZOSH];
    for (int i=0; i<nOriginalCount; i++)
    {
        float fValue = fCrossResult[i];
        int iIndex = i;
        for (int j=i + 1; j<nOriginalCount; j++)
        {
            if (fCrossResult[j] < fValue)
            {
                fValue = fCrossResult[j];
                iIndex = j;
            }
        }

        if (iIndex != i)
        {
            std::swap(fCrossResult[iIndex], fCrossResult[i]);
            std::swap(iCrossIndices[iIndex], iCrossIndices[i]);
        }
    }


    //   说明：预处理索引，极角相同时取最远的 [10/21/2016 ZOSH];
    int iHandleIndices[MAX_POINT_NUMBER] = {0};
    int nHandleCount = 0;
    //   说明：先将凸点加入 [10/21/2016 ZOSH];
    iHandleIndices[nHandleCount++] = iStartPointIndex;
    for (int i=0; i<nOriginalCount; i++)
    {
        if (iCrossIndices[i] == iStartPointIndex)
        {
            continue;
        }

        float fValue_i = fCrossResult[i];
        iHandleIndices[nHandleCount] = iCrossIndices[i];

        for (int j=i+1; j<nOriginalCount; j++)
        {
            float fValue_j = fCrossResult[j];

            if ( fabs(fValue_j - fValue_i) < 1e-6)
            {
                //int iIndex_i = iCrossIndices[i];
                int iIndex_i = iHandleIndices[nHandleCount];
                int iIndex_j = iCrossIndices[j];

                //   说明：计算距离 [10/21/2016 ZOSH];
                int iSqDist_i = int(pow(pOriginalPos[iIndex_i].x() - pOriginalPos[iStartPointIndex].x(), 2) + pow(pOriginalPos[iIndex_i].y() - pOriginalPos[iStartPointIndex].y(), 2));
                int iSqDist_j = int(pow(pOriginalPos[iIndex_j].x() - pOriginalPos[iStartPointIndex].x(), 2) + pow(pOriginalPos[iIndex_j].y() - pOriginalPos[iStartPointIndex].y(), 2));

                if (iSqDist_j > iSqDist_i)
                {
                    iHandleIndices[nHandleCount] = iIndex_j;
                }

                continue;
            }
            else
            {
                int iStepCount = j;
                i+= iStepCount - i - 1;
                break;
            }
        }


        nHandleCount++;
    }

    //   说明：个数小于2， 直接返回 [10/21/2016 ZOSH];
    if (nHandleCount < 2)
    {
        return -1;
    }

    std::vector<int> stackConvexIndices;
    for (int i=0; i<3; i++)
    {
        stackConvexIndices.push_back(iHandleIndices[i]);
    }

    for (int i=3; i<nHandleCount; i++)
    {
        int iCurrentIndex = iHandleIndices[i];

        while (true)
        {
            if (stackConvexIndices.size() > 1)
            {
                int iTopIndex = stackConvexIndices.back();
                stackConvexIndices.pop_back();
                int iNextTopIndex = stackConvexIndices.back();

                int i_value_pi_ptop[2] = {pOriginalPos[iTopIndex].x() - pOriginalPos[iNextTopIndex].x(), pOriginalPos[iTopIndex].y() - pOriginalPos[iNextTopIndex].y()};
                int i_value_pi_pnext_top[2] = { pOriginalPos[iCurrentIndex].x() - pOriginalPos[iNextTopIndex].x(), pOriginalPos[iCurrentIndex].y() - pOriginalPos[iNextTopIndex].y()};
                //QVector3D i_value_pi_ptop(pOriginalPos[iTopIndex].x() - pOriginalPos[iNextTopIndex].x(), pOriginalPos[iTopIndex].y() - pOriginalPos[iNextTopIndex].y(), 0);
                //QVector3D i_value_pi_pnext_top(pOriginalPos[iCurrentIndex].x() - pOriginalPos[iNextTopIndex].x(), pOriginalPos[iCurrentIndex].y() - pOriginalPos[iNextTopIndex].y(), 0);

                //int iValueCrossDot = QVector3D::crossProduct(i_value_pi_ptop,i_value_pi_pnext_top);
                int iValueCrossDot = i_value_pi_ptop[0]*i_value_pi_pnext_top[1] - i_value_pi_ptop[1]*i_value_pi_pnext_top[0];

                //   说明：注意这里方向可能错误 [10/21/2016 ZOSH];
                if (iValueCrossDot > 0)
                {
                    stackConvexIndices.push_back(iTopIndex);
                    stackConvexIndices.push_back(iCurrentIndex);
                    break;
                }
                else
                {
                    continue;
                };
            }
            else
            {
                break;
            }
        }
    }

    while (!stackConvexIndices.empty())
    {
        iConvexIndices[nConvexPointCount++] = stackConvexIndices.back();
        stackConvexIndices.pop_back();
    }

    return 0;
}

void BorderFinder::creatConvexHull(string intputFile, string outputFile)
{
    GrahamScan ex;
    vector<Node> stack = ex.readFile(intputFile);
    Node nodeSaved = ex.outputFile(outputFile);
    QFileInfo fileInfo(QString::fromStdString(outputFile));
    QString fileName = fileInfo.baseName();
    QString filepng(QCoreApplication::applicationDirPath() + + "/Borderfinder/" + fileName + ".png");
    QVector<QPointF> pointList;
    QPointF point;
    for(auto node : stack)
    {
        point.setX(node.getX()+nodeSaved.getX());
        point.setY(node.getY()+nodeSaved.getY());
        pointList.push_back(point);
    }
    qDebug()<<"creatConvexHull : "<<pointList.size();
    qDebug()<<"filepng : "<<filepng;
    qDebug()<<"pointList : "<<pointList;
    drawBorderPointsToImage(pointList,filepng, true);
}

bool BorderFinder::collisionDetection(string filename1, string filename2)
{
    ifstream myFile;
    string temp;
    int numNodes, lineNum = 1;

    QTextCodec *code = QTextCodec::codecForName("GBK");
    string str = code->fromUnicode(QString::fromStdString(filename1)).data();
    //cout << "str..." << str << endl;
    myFile.open(str.c_str());
    if (!myFile)
    {
        cout << "File Failed to Open! Filename:" << str << endl;
        exit(1);
    }
    cout << "Reading file..." << endl;
    //Get the file from input
    getline(myFile, temp);
    numNodes = atoi(temp.c_str()); //The first line is always the number of Nodes in the graph
    vector<Node> graph1;
    while (getline(myFile, temp))
    {
        int comma = temp.find(",");
        string xVal = temp.substr(0, comma);
        string yVal = temp.substr(comma + 1);
        graph1.push_back(Node(atoi(xVal.c_str()), atoi(yVal.c_str()))); //adding the new node to the graph
    }
    myFile.close();
    //qDebug()<<"graph1: "<<graph1.size();

    str = code->fromUnicode(QString::fromStdString(filename2)).data();
    //cout << "str..." << str << endl;
    myFile.open(str.c_str());
    if (!myFile)
    {
        cout << "File Failed to Open! Filename:" << str << endl;
        exit(1);
    }
    cout << "Reading file..." << endl;
    //Get the file from input
    getline(myFile, temp);
    numNodes = atoi(temp.c_str()); //The first line is always the number of Nodes in the graph
    vector<Node> graph2;
    while (getline(myFile, temp))
    {
        int comma = temp.find(",");
        string xVal = temp.substr(0, comma);
        string yVal = temp.substr(comma + 1);
        graph2.push_back(Node(atoi(xVal.c_str()), atoi(yVal.c_str()))); //adding the new node to the graph
    }
    myFile.close();
    //qDebug()<<"graph2: "<<graph2.size();

    bool bDetection = false;
    for(int i =0;i<graph1.size();i++)
    {
        for(int j =0;j<graph2.size();j++)
        {
            if ( max(graph1[i].getX(), graph1[(i+1)%(graph1.size())].getX())<min(graph2[j].getX(), graph2[(j+1)%(graph2.size())].getX()) )
            {
                //qDebug()<<"collisionDetection1: "<< "max: " << max(graph1[i].getX(), graph1[(i+1)%(graph1.size())].getX()) << ", "<<"min: " << min(graph2[j].getX(), graph2[(j+1)%(graph2.size())].getX());
                bDetection = false;
            }
            else if ( max(graph1[i].getY(), graph1[(i+1)%(graph1.size())].getY())<min(graph2[j].getY(), graph2[(j+1)%(graph2.size())].getY()) )
            {
                //qDebug()<<"collisionDetection2";
                bDetection = false;
            }
            else if ( max(graph2[i].getX(), graph2[(i+1)%(graph2.size())].getX())<min(graph1[j].getX(), graph1[(j+1)%(graph1.size())].getX()) )
            {
                //qDebug()<<"collisionDetection3: "<< "max: " << max(graph2[i].getX(), graph2[(i+1)%(graph2.size())].getX()) << ", "<<"min: " << min(graph1[j].getX(), graph1[(j+1)%(graph1.size())].getX());
                bDetection = false;
            }
            else if ( max(graph2[i].getY(), graph2[(i+1)%(graph2.size())].getY())<min(graph1[j].getY(), graph1[(j+1)%(graph1.size())].getY()) )
            {
                //qDebug()<<"collisionDetection4";
                bDetection = false;
            }
            else if ( mult(graph2[j], graph1[(i+1)%(graph1.size())], graph1[i])*mult(graph1[(i+1)%(graph1.size())], graph2[(j+1)%(graph2.size())], graph1[i])<0 )
            {
                //qDebug()<<"collisionDetection5";
                bDetection = false;
            }
            else if ( mult(graph1[i], graph2[(j+1)%(graph2.size())], graph2[j])*mult(graph2[(j+1)%(graph2.size())], graph1[(i+1)%(graph1.size())], graph2[j])<0 )
            {
                //qDebug()<<"collisionDetection6";
                bDetection = false;
            }
            else
            {
                //qDebug()<<"collisionDetection7";
                return true;
            }
        }
    }
    return bDetection;
}

double BorderFinder::mult(Node a, Node b, Node c)
{
    return (a.getX()-c.getX())*(b.getY()-c.getY())-(b.getX()-c.getX())*(a.getY()-c.getY());
}

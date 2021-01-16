#include <QCoreApplication>
#include <qdir.h>
#include "borderfinder.h"
#include <QTextStream>
#include <QDebug>
#include <QPainter>
#include <QTextCodec>
#include "borderfinderConstants.h"
#include <QTextCodec>


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
        qDebug()<<"read_mesh";
        if (!OpenMesh::IO::read_mesh(myMesh, str))
        {
            qDebug()<<"src_path::"<<src_path;
            std::cout << "Fail to read mesh!" << std::endl;
            return false;
        }
        myMeshBuff = myMesh;
    }
    //QString tempPath1 = QCoreApplication::applicationDirPath()+"/Borderfinder/C1.stl";
    string tempPath1 = code->fromUnicode(QCoreApplication::applicationDirPath()+"/Borderfinder/C1.stl").data();
    double time_Start = (double)clock();
    OpenMesh::IO::write_mesh(myMesh, tempPath1, OpenMesh::IO::Options::Binary);
    OpenMesh::IO::read_mesh(myMesh, tempPath1);
    double time_End = (double)clock();
    qDebug()<<"time: "<<(time_End - time_Start)/1000.0<<"s";

    qDebug()<<"111";
    //myMesh.halfedges();
    //myMesh.set_halfedge_handle();
    myMesh.request_vertex_status();
    myMesh.request_halfedge_status();
    myMesh.request_edge_status();
    myMesh.request_face_status();
    myMesh.request_face_normals();
    myMesh.update_normals();
    const MyPoint standard_normal(0, 0, 1);
    qDebug()<<"222";
    int i = 0;
    int k = 0;
    for (auto f_it = myMesh.faces_begin(); f_it != myMesh.faces_end(); f_it++)
    {
        k++;
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
        //qDebug()<<"x = "<<x<<","<<"y = "<<y<<","<<"z = "<<z;
        if ((x*standard_normal[0] + y * standard_normal[1] + z * standard_normal[2]) < 0)
        {
                    i++;
            myMesh.delete_face(*f_it);
        }
    }
    qDebug()<<"myMesh.facce = "<<i;
    qDebug()<<"myMesh.faccek = "<<k;
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
    qDebug()<<"myMesh.halfedges = "<<j;
    qDebug()<<"444";
    ClipperLib::Clipper clipper;
    ClipperLib::Paths src_paths;
    ClipperLib::Paths dst_paths;
    ClipperLib::ClipperOffset  clipperoff;
    const double scale = std::pow(10, precision);
    qDebug()<<"myMesh.boundaries = "<<boundaries.size();
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
//    qDebug()<<"src_paths.size = "<<src_paths.size();
//    qDebug()<<"444444";
    //ClipperLib::SimplifyPolygons(src_paths,src_paths,ClipperLib::pftNonZero);
    clipper.AddPaths(src_paths, ClipperLib::PolyType::ptSubject, true);
//    qDebug()<<"4443444";
    clipper.Execute(ClipperLib::ClipType::ctUnion, dst_paths, ClipperLib::PolyFillType::pftNonZero);
//    qDebug()<<"44433444";
    ClipperLib::CleanPolygons(dst_paths,10);
//    qDebug()<<"444333444";
    clipperoff.AddPaths(dst_paths,ClipperLib::jtSquare,ClipperLib::etClosedPolygon);
//    qDebug()<<"4443333444";
    clipperoff.Execute(dst_paths,0 * scale);
//    qDebug()<<"44433333444";
    ClipperLib::SimplifyPolygons(dst_paths,dst_paths, ClipperLib::PolyFillType::pftNonZero);
//    qDebug()<<"444333333444";
    boundaries.clear();
//    qDebug()<<"555";
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

    qDebug()<<"boundaries.size = : "<<boundaries.size();
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
    qDebug()<<"renderScale = : "<<renderScale;
    for (const auto& boundary : boundaries)
    {
        for (const auto& point : boundary)
        {           
            QPointF tempPoint(point[0]*(1.0/renderScale),point[1]*(1.0/renderScale));
            pointList1.push_back(tempPoint);
        }
        //break;
    }
    if(!drawBorderPointsToImage(pointList1,filepng))
    {
        qDebug()<<"drawBorderPointsToImage: false";
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

        qDebug()<<"pointList.size()"<<pointList.size();

        out << pointList.size();
        for (QPointF point : pointList) {
             out << ' ' << '[' << point.x() << ',' << 2*coorH - point.y() << ']';
            //out << ' ' << '[' << point.x() << ',' << point.y() + y_offset<< ']';
        }
        out << endl;
        file.close();
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
    qDebug()<<"cloneModelAndProcess1";
    qDebug()<<"cloneModelAndProcess1: "<<QString::fromStdString(inputPath);
    QTextCodec *code = QTextCodec::codecForName("GBK");
    string str = code->fromUnicode(QString::fromStdString(inputPath)).data();
    if (!OpenMesh::IO::read_mesh(myMesh, str))
    {
        std::cout << str << std::endl;
        qDebug()<<"Failed to read mesh!";
        std::cout << "Failed to read mesh!" << std::endl;
        return false;
    }

    qDebug()<<"cloneModelAndProcess2";
    offX = 0, offY = 0;
    bool isNormal = calcOff(offX, offY, scale);
    if(!isNormal)
    {
        return isNormal;
    }

    QString generatePath1 = QCoreApplication::applicationDirPath() + BORDERFINDERCONSTANTS::CLONE_STL_PATH + "/d.stl";
    auto localPath1 = generatePath1.toLocal8Bit();
    qDebug()<<"read_mesh11111: "<<QString::fromStdString(inputPath);
    qDebug()<<"write_mesh11111: "<<localPath1;
    str = code->fromUnicode(localPath1).data();
    OpenMesh::IO::write_mesh(myMesh, str, OpenMesh::IO::Options::Binary);

    qDebug()<<"cloneModelAndProcess3";
    if(offX != 0 || offY != 0)
    {
        QTextCodec *code = QTextCodec::codecForName("GBK");
        string str = code->fromUnicode(QString::fromStdString(inputPath)).data();
        QString qinputPath = QString(QString::fromLocal8Bit(str.c_str()));
        qDebug()<<"qinputPath: "<<qinputPath;
        QDir dir(QCoreApplication::applicationDirPath() + BORDERFINDERCONSTANTS::SQUASHED_IMG_PATH);
        QFileInfo fileInfo(qinputPath);
        QString fileName = fileInfo.baseName();
        //QString fileName = qinputPath.right(qinputPath.size()-qinputPath.lastIndexOf('.')-1);
        qDebug()<<"qinputPath.fileName: "<<fileName;
        qDebug()<<"modelFilePath: "<<modelFilePath;
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
        qDebug()<<"generatePath: "<<generatePath;
        //auto localPath = generatePath.toLocal8Bit();
        //qDebug()<<"read_mesh11111: "<<QString::fromStdString(inputPath);
        //qDebug()<<"write_mesh11111: "<<localPath;
        str = code->fromUnicode(generatePath).data();
        if (!OpenMesh::IO::write_mesh(myMesh, str, OpenMesh::IO::Options::Binary))
        {
            std::cerr << "Cannot write mesh to file ' .stl ' " << std::endl;
            qDebug()<<"write_mesh err";
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
        qDebug()<<"myMesh.vertices_empty";
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

bool BorderFinder::drawBorderPointsToImage(QVector<QPointF> pointList,QString imagePath)
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

//    double offsetX = 0;
//    double offsetY = 0;
    double y_offset = 1080 - 2 * coorH;

//    if(minX < 0 || minY < 0)
//    {
//       return false;
//    }

    ramerDouglasPeucker(clonePointList, 15.0, pointListOut);

    QPainter p(&borderedImg);
    p.setPen(QColor(255, 255, 255));
    int length = pointListOut.size();
    qDebug()<<"pointListOut.size: "<<length;
    for (int i = 0; i < length; i++)
    {
        //qDebug()<<"(x1,y1)= "<<"(" <<pointListOut[i%length].x() - offsetX<<","<<pointListOut[i%length].y() - offsetY<<")"<<endl;
        //qDebug()<<"(x2,y2)= "<<"(" <<pointListOut[(i + 1) % length].x() - offsetX<<","<<pointListOut[(i + 1) % length].y() - offsetY<<")"<<endl;
        //p.drawLine(pointListOut[i%length].x() - offsetX, pointListOut[i%length].y() - offsetY, pointListOut[(i + 1) % length].x() - offsetX, pointListOut[(i + 1) % length].y() - offsetY);
        p.drawLine(pointListOut[i%length].x(), 2*coorH - pointListOut[i%length].y(), pointListOut[(i + 1) % length].x(), 2*coorH - pointListOut[(i + 1) % length].y());
        //p.drawLine(pointListOut[i%length].x(), pointListOut[i%length].y()+y_offset, pointListOut[(i + 1) % length].x(), pointListOut[(i + 1) % length].y()+y_offset);
        //p.drawLine(pointListOut[i%length].x(), pointListOut[i%length].y(), pointListOut[(i + 1) % length].x(), pointListOut[(i + 1) % length].y());
    }

    // will return false if error(point.y() + y_offset)
    return borderedImg.save(imagePath);
}

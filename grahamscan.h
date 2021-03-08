#ifndef GRAHAMSCAN_H
#define GRAHAMSCAN_H

#include <vector>
/**modified by wangjx Node.h->node.h
  *it will find "node.h" to match Node.h,but that not work on mac .So change it.
  */
#include "node.h"

using namespace std;
class GrahamScan
{
private:
    vector<Node> stack; //For hull-making
    vector<Node> graph;
    int numOfNodes;
    int savedX, savedY;

public:
    GrahamScan();
    ~GrahamScan();
    vector<Node> readFile(string filename); //reads the file and calls required execution
    Node outputFile(string filename); //outputs file named   original_hull.txt
    void turnDirection(Node prev1st, Node prev2nd, Node prev3rd); //Calculates which turn we are taking
    void setNumOfNodes(int num);
    int getNumOfNodes();
    void scan(); //This is the main event
    void findLowestY();
};

#endif // GRAHAMSCAN_H

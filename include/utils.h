#ifndef UTILS_H
#define UTILS_H

#include "point3d.h"
#include <vector>



std::vector<Point3d> generate3dPoints(int n, double xmin, double xmax, double ymin, double ymax,  double zmin, double zmax);
/*

  {

  std::vector<double> vecOfRandomNumsX(n);
  std::vector<double> vecOfRandomNumsY(n);
  std::vector<double> vecOfRandomNumsZ(n);
  
  srand(0);
  
  std::vector<Point3d> points;
  
  std::generate(vecOfRandomNumsX.begin(), vecOfRandomNumsX.end(), []()
  {
    return static_cast <double> (rand() % 2001)/2000.0;
  });
  
  std::generate(vecOfRandomNumsY.begin(), vecOfRandomNumsY.end(), []()
  {
    return static_cast <double> (rand() % 2001)/2000.0;
  });
  
  std::generate(vecOfRandomNumsZ.begin(), vecOfRandomNumsZ.end(), []()
  {
    return static_cast <double> (rand() % 2001)/2000.0;
  });
  
  for (int i = 0; i < n; i++) {
    points[i].x = vecOfRandomNumsX[i]*(xmax-xmin) + xmin;
    points[i].y = vecOfRandomNumsY[i]*(ymax-ymin) + ymin;
    points[i].z = vecOfRandomNumsZ[i]*(zmax-zmin) + zmin;
  }
  
  return points;
  

  
}
*/


#endif

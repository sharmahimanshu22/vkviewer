#ifndef UTILS_H
#define UTILS_H


#include "point3d.h"
#include <vector>


inline std::vector<Point3d> generate3dPoints(int n, double xmin, double xmax, double ymin, double ymax,  double zmin, double zmax, int seed)

  {

  std::vector<double> vecOfRandomNumsX(n);
  std::vector<double> vecOfRandomNumsY(n);
  std::vector<double> vecOfRandomNumsZ(n);
  
  srand(0);
  
  std::vector<Point3d> points(n);
  
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

/*
template <class T>
struct cmp_zsort_firstpoint {
public:
  cmp_zsort_firstpoint(uint8_t idx) : idx_(idx) {}
  
  bool operator() (const T& a, const T& b) const
  {
    // Custom comparison logic
    return vertices[a[0]][idx_] < vertices[b[0]][idx_]; // this sorts in ascending order
  }
  
private:
  int idx_;
};
  


template <class RandomIt, class Compare>
inline void zsort_triangles(RandomIt first, RandomIt last, Compare comp) {

  std::queue<std::tuple<RandomIt, RandomIt, uint8_t>> sort_tuples;
  sort_range.push(std::make_tuple(first, last, 0));

  cmp_zsort cmpx(0);
  cmp_zsort cmpy(1);
  cmp_zsort cmpz(2);
  cmp_zsort cmps[3] = {cmpx, cmpy, cmpz};
  
  while(!sort_tuples.empty()) {
    std::tuple<RandomIt, RandomIt, uint8_t> element = sort_tuples.front();
    sort_tuples.pop();
    RandomIt first = e[0];
    RandomIt last = e[1];
    if(first == last) continue;
    uint8_t idx = e[2];

    std::size = std::distance(first, last)/2;
    std::nth_element(first, first + mid , last, cmps[idx]);

    uint8_t nidx = (idx + 1)%3;
    sort_tuples.push(std::make_tuple(first, mid, nidx));
    sort_tuples.push(std::make_tuple(mid+1, last, nidx));
  }
  
  // it should be sorted now

}
*/

#endif

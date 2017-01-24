#include "CollisionDetection.h"
#include "CrispCollisionDetector.h"
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <itksys/SystemTools.hxx>
#include "itkExceptionObject.h"

int main(int argc, char** argv)
{
  CrispCollisionDetector detector(argc, argv);
  return 0;
}

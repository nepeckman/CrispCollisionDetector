#ifndef CRISP_COLLISION_DETECTOR_H
#define CRISP_COLLISION_DETECTOR_H
#include <Eigen/Dense>
#include <unordered_set>
#include "CollisionDetection.h"
#include "ItkTypes.h"



struct Link;
struct Node{
	Image3DType::IndexType *index;
	Link *links;
	Node *next;
	Node *previous;
	int processed;
};

struct Link{
	Node *to;
	Node *from;
	Link *next;
};

struct Line3D {
	Eigen::Vector3d pt1;
	Eigen::Vector3d pt2;
};

struct Cylinder{
	Line3D* line;
	float radius;
};

struct Cube{
  Eigen::Vector3d center;
  float halfSide;
};

typedef std::vector< std::unordered_set< Node* >* > ObjectVectorType;

class CrispCollisionDetector : public CollisionDetection {
public:
	Image3DType* image;
	ReaderType::Pointer reader;
	ImageIOType::Pointer dicomIO;
	Image3DType::SizeType size;
	bool inCollision(const std::vector<Eigen::Vector3d> & point1s,
			   					 const std::vector<Eigen::Vector3d> & point2s,
			   			 		 const std::vector<float> & radii,
			   			 		 std::vector<int> & indices,
								 	 std::vector<Image3DType::IndexType*> & pixels) const;
	void remove_starting_points(const std::vector<Eigen::Vector3d> & point1s,
															const std::vector<Eigen::Vector3d> & point2s,
															const std::vector<float> & radii);
	CrispCollisionDetector(int argc, char** argv);

	bool ptIsInCylinder(const Eigen::Vector3d &pt, Cylinder* cylinder) const;
	void visualize(int argc, char** argv);
	Eigen::Vector3d ptOnLine(const Eigen::Vector3d &pt, Line3D* line) const;
	Eigen::Vector3d moveAlongLine(Line3D* line, float distance) const;
	float ptToLineDistance(const Eigen::Vector3d &pt, Line3D* line) const;
	bool pointIsOnSegment(const Eigen::Vector3d &pt, Line3D *line) const;
	float magnitude(const Eigen::Vector3d &pt) const;
	float pointDistance(const Eigen::Vector3d &pt1, const Eigen::Vector3d &pt2) const;
	bool distanceEqual(float d1, float d2) const;
	int findIndex(Image3DType::SizeType size, Image3DType::IndexType index) const;
	int findIndex(Image3DType::SizeType size, int x, int y, int z) const;
	int isEqual(Image3DType::IndexType index1, Image3DType::IndexType index2) const;
	float cylinderMaxDistance(Cylinder* cylinder) const;
	Cube* boundingCube(Cylinder* cylinder) const;
	Eigen::Vector3d lowestPoint(Cube* cube) const;
	Eigen::Vector3d highestPoint(Cube* cube) const;
  std::vector< Node* > getGraph();
};
#endif

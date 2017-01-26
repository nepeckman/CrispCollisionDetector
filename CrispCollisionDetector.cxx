#include "CrispCollisionDetector.h"
#include "CollisionDetection.h"
#include "ItkTypes.h"
#include <iostream>
#include <string>
#include <vector>
#include <unordered_set>
#include <Eigen/Dense>
#include <queue>
#include <algorithm>
#include <cmath>
#include <itksys/SystemTools.hxx>
#include "itkExceptionObject.h"
#include <limits>

// By Thursday:
// Switch doubles to floats
// Remove configurable cube from starting points
// Testing
CrispCollisionDetector::CrispCollisionDetector(int argc, char** argv)
{
	// create name generator
	std::cout << "Creating name generator" << std::endl;
	NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
	nameGenerator->SetDirectory("../../PEL");

	// create dicom reader
	std::cout << "Creating dicom reader" << std::endl;
	reader = ReaderType::New();
	dicomIO = ImageIOType::New();
	reader->SetImageIO(dicomIO);

	// get series IDs
	std::cout << "Creating dicom series ID" << std::endl;
	const SeriesIDContainer &seriesUID = nameGenerator->GetSeriesUIDs();

	// get dicom series
	std::cout << "Getting dicom series" << std::endl;
	FileNamesContainer fileNames;
	fileNames = nameGenerator->GetFileNames(seriesUID.begin()->c_str());
	reader->SetFileNames(fileNames);
	reader->Update();

	// read image data
	std::cout << "Reading image data" << std::endl;
	image = reader->GetOutput();
	Image3DType::RegionType region = image->GetLargestPossibleRegion();
	size = region.GetSize();

	std::cout << "Finding max indecies" << std::endl;
	int maxValue = std::numeric_limits<int>::max();
	int minIndicies[3] = {maxValue, maxValue, maxValue};
	int maxIndicies[3] = {0, 0, 0};
	for(int i = 0; i < size[0]; i++)
	{
		for(int j = 0; j < size[1]; j++)
		{
			for(int k = 0; k < size[2]; k++)
			{
				Image3DType::IndexType index;
				index[0] = i;
				index[1] = j;
				index[2] = k;
				Image3DType::PixelType pixelValue = image->GetPixel(index);
				if(pixelValue > 0)
				{
					image->SetPixel(index, 0);
					// find the bounding cube
					for(int idx = 0; idx < 3; idx++)
					{
						if(index[idx] > maxIndicies[idx])
						{
							maxIndicies[idx] = index[idx];
						}
						if(index[idx] < minIndicies[idx])
						{
							minIndicies[idx] = index[idx];
						}
					}
				}
			}
		}
	}

	std::cout << "Filling in binding cube" << std::endl;
	int wall = 3;
	for(int i = minIndicies[0] - wall; i <= minIndicies[0]; i++ )
	{
			for(int j = minIndicies[1] - wall; j <= maxIndicies[1] + wall; j++)
			{
				for(int k = minIndicies[2] - wall; k <= maxIndicies[2] + wall; k++)
				{
					Image3DType::IndexType index;
					index[0] = i;
					index[1] = j;
					index[2] = k;
					image->SetPixel(index, 10000);
				}
			}
	}

for(int i = maxIndicies[0]; i <= maxIndicies[0] + wall; i++ )
	{
			for(int j = minIndicies[1] - wall; j <= maxIndicies[1] + wall; j++)
			{
				for(int k = minIndicies[2] - wall; k <= maxIndicies[2] + wall; k++)
				{
					Image3DType::IndexType index;
					index[0] = i;
					index[1] = j;
					index[2] = k;
					image->SetPixel(index, 10000);
				}
			}
	}

for(int i = minIndicies[0] - wall; i <= maxIndicies[0] + wall; i++ )
	{
			for(int j = minIndicies[1] - wall; j <= minIndicies[1]; j++)
			{
				for(int k = minIndicies[2] - wall; k <= maxIndicies[2] + wall; k++)
				{
					Image3DType::IndexType index;
					index[0] = i;
					index[1] = j;
					index[2] = k;
					image->SetPixel(index, 10000);
				}
			}
	}

for(int i = minIndicies[0] - wall; i <= maxIndicies[0] + wall; i++ )
	{
			for(int j = maxIndicies[1]; j <= maxIndicies[1] + wall; j++)
			{
				for(int k = minIndicies[2] - wall; k <= maxIndicies[2] + wall; k++)
				{
					Image3DType::IndexType index;
					index[0] = i;
					index[1] = j;
					index[2] = k;
					image->SetPixel(index, 10000);
				}
			}
	}

for(int i = minIndicies[0] - wall; i <= maxIndicies[0] + wall; i++ )
	{
			for(int j = minIndicies[1] - wall; j <= maxIndicies[1] + wall; j++)
			{
				for(int k = minIndicies[2] - wall; k <= minIndicies[2]; k++)
				{
					Image3DType::IndexType index;
					index[0] = i;
					index[1] = j;
					index[2] = k;
					image->SetPixel(index, 10000);
				}
			}
	}

for(int i = minIndicies[0] - wall; i <= maxIndicies[0] + wall; i++ )
	{
			for(int j = minIndicies[1] - wall; j <= maxIndicies[1] + wall; j++)
			{
				for(int k = maxIndicies[2]; k <= maxIndicies[2] + wall; k++)
				{
					Image3DType::IndexType index;
					index[0] = i;
					index[1] = j;
					index[2] = k;
					image->SetPixel(index, 10000);
				}
			}
	}

	// make output directory
	std::string writedir;
		if(argc < 2)
		{
			writedir = "out";
		} else
		{
			writedir = argv[2];
		}
		std::cout << "Write directory is \"" << writedir << "\"" << std::endl;


		itksys::SystemTools::MakeDirectory(writedir);

		// generate the file names
		OutputNamesGeneratorType::Pointer outputNames = OutputNamesGeneratorType::New();
		std::string seriesFormat = std::string(writedir) + "/image-%05d.dcm";
		outputNames->SetSeriesFormat(seriesFormat.c_str());
		outputNames->SetStartIndex(1);
		outputNames->SetEndIndex(size[2]);

		// create series writer
		reader->Update();

		SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
		seriesWriter->SetInput(image);
		seriesWriter->SetImageIO(dicomIO);
		seriesWriter->SetFileNames(outputNames->GetFileNames());
		seriesWriter->SetMetaDataDictionaryArray(reader->GetMetaDataDictionaryArray());

		try
		{
			seriesWriter->Update();
		}
		catch(itk::ExceptionObject & excp)
		{
			std::cerr << "Exception throw while writing the series" << std::endl;
			std::cerr << excp << std::endl;
		}
}

int CrispCollisionDetector::isEqual(Image3DType::IndexType index1, Image3DType::IndexType index2) const
{
	return (index1[0] == index2[0] && index1[1] == index2[1] && index1[2] == index2[2]);
}

int CrispCollisionDetector::findIndex(Image3DType::SizeType size, Image3DType::IndexType index) const
{
	return index[0] * size[1] * size[2] + index[1] * size[2] + index[2];
}

int CrispCollisionDetector::findIndex(Image3DType::SizeType size, int x, int y, int z) const
{
	return x * size[1] * size[2] + y * size[2] + z;
}

bool CrispCollisionDetector::distanceEqual(float d1, float d2) const
{
	return std::abs(d1 - d2) < 0.01;
}

float CrispCollisionDetector::pointDistance(const Eigen::Vector3d &pt1, const Eigen::Vector3d &pt2) const
{
	float diffx = std::abs(pt1(0) - pt2(0));
	float diffy = std::abs(pt1(1) - pt2(1));
	float diffz = std::abs(pt1(2) - pt2(2));
	return std::sqrt(std::pow(diffx, 2) + std::pow(diffy, 2) + std::pow(diffz, 2));
}

float CrispCollisionDetector::magnitude(const Eigen::Vector3d &pt) const
{
	return std::sqrt(std::pow(pt(0), 2) + std::pow(pt(1), 2) + std::pow(pt(2), 2));
}

bool CrispCollisionDetector::pointIsOnSegment(const Eigen::Vector3d &pt, Line3D *line) const
{
	return distanceEqual(pointDistance(pt, line->pt1) + pointDistance(pt, line->pt2), pointDistance(line->pt1, line->pt2));
}

float CrispCollisionDetector::ptToLineDistance(const Eigen::Vector3d &pt, Line3D* line) const
{
	float numerator = magnitude((pt - line->pt1).cross(pt - line->pt2));
	float denominator = magnitude(line->pt2 - line->pt1);
	return numerator/denominator;
}

Eigen::Vector3d CrispCollisionDetector::moveAlongLine(Line3D* line, float distance) const
{
	Eigen::Vector3d vector = line->pt2 - line->pt1;
	Eigen::Vector3d movedVector = vector * distance/pointDistance(line->pt1, line->pt2);
	return line->pt1 + movedVector;
}

Eigen::Vector3d CrispCollisionDetector::ptOnLine(const Eigen::Vector3d &pt, Line3D* line) const
{
	float ptToSeg = ptToLineDistance(pt, line);
	float linePtToIntersection = std::sqrt(std::pow(pointDistance(line->pt1, pt), 2) - std::pow(ptToSeg, 2));
	return moveAlongLine(line, linePtToIntersection);
}

bool CrispCollisionDetector::ptIsInCylinder(const Eigen::Vector3d &pt, Cylinder* cylinder) const
{
	Eigen::Vector3d intersection = ptOnLine(pt, cylinder->line);
	if(pointIsOnSegment(intersection, cylinder->line))
	{
		return pointDistance(intersection, pt) < cylinder->radius;
	} else
	{
		return false;
	}
}

float CrispCollisionDetector::cylinderMaxDistance(Cylinder* cylinder) const
{
	Eigen::Vector3d midpoint = (cylinder->line->pt1 + cylinder->line->pt2)/2;
	float a = pointDistance(midpoint, cylinder->line->pt2);
	float b = cylinder->radius;
	return std::sqrt(std::pow(a, 2) + std::pow(b, 2));
}

Cube* CrispCollisionDetector::boundingCube(Cylinder* cylinder) const
{
	float halfSide = cylinderMaxDistance(cylinder);
	Eigen::Vector3d midpoint = (cylinder->line->pt1 + cylinder->line->pt2)/2;
	Cube* box = new Cube;
	box->halfSide = halfSide;
	box->center = midpoint;
	return box;
}

Eigen::Vector3d CrispCollisionDetector::lowestPoint(Cube* cube) const
{
	Eigen::Vector3d lowest;
	lowest(0) = cube->center(0) - cube->halfSide;
	lowest(1) = cube->center(1) - cube->halfSide;
	lowest(2) = cube->center(2) - cube->halfSide;
	return lowest;
}
Eigen::Vector3d CrispCollisionDetector::highestPoint(Cube* cube) const
{
	Eigen::Vector3d highest;
	highest(0) = cube->center(0) + cube->halfSide;
	highest(1) = cube->center(1) + cube->halfSide;
	highest(2) = cube->center(2) + cube->halfSide;
	return highest;
}

bool CrispCollisionDetector::inCollision(const std::vector<Eigen::Vector3d> & point1s,
								 const std::vector<Eigen::Vector3d> & point2s,
								 const std::vector<float> & radii,
								 std::vector<int> & indices,
							 	 std::vector<Image3DType::IndexType*> & pixels) const
{
		bool rv = false;
		for(std::vector< Eigen::Vector3d >::size_type i = 0; i != point1s.size(); i++)
		{
			Cylinder* cylinder = new Cylinder;
			cylinder->line = new Line3D;
			cylinder->line->pt1 = point1s[i];
			cylinder->line->pt2 = point2s[i];
			cylinder->radius = radii[i];
			Cube* cube = boundingCube(cylinder);
			Eigen::Vector3d lowest = lowestPoint(cube);
			Eigen::Vector3d highest = highestPoint(cube);
			int x1 = (int) std::floor(lowest(0));
			int y1 = (int) std::floor(lowest(1));
			int z1 = (int) std::floor(lowest(2));
			int x2 = (int) std::ceil(highest(0));
			int y2 = (int) std::ceil(highest(1));
			int z2 = (int) std::ceil(highest(2));
			for(int x = x1; x <= x2; x++)
			{
				for(int y = y1; y <= y2; y++)
				{
					for(int z = z1; z <= z2; z++)
					{
						Eigen::Vector3d point;
						point(0) = x;
						point(1) = y;
						point(2) = z;
						Image3DType::IndexType *pixel = new Image3DType::IndexType;
						(*pixel)[0] = x;
						(*pixel)[1] = y;
						(*pixel)[2] = z;
						if(ptIsInCylinder(point, cylinder) && (image->GetPixel(*pixel) > 0))
						{
							rv = true;
							indices.push_back(i);
							pixels.push_back(pixel);
						}
					}
				}
			}
			delete cylinder->line;
			delete cylinder;
		}
		return false;
}

void CrispCollisionDetector::remove_starting_points(const std::vector<Eigen::Vector3d> & point1s,
								 const std::vector<Eigen::Vector3d> & point2s,
								 const std::vector<float> & radii)
{
	std::vector<int> indices;
	std::vector<Image3DType::IndexType*> pixels;
	inCollision(point1s, point2s, radii, indices, pixels);
	for(std::vector<Image3DType::IndexType*>::size_type i = 0; i != pixels.size(); i++)
	{
		image->SetPixel(*pixels[i], 0);
	}
	return;
}

#ifndef ITK_TYPES_H
#define ITK_TYPES_H
#include "itkImage.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkThresholdImageFilter.h"
#include "itkImageFileWriter.h"
#include <itksys/SystemTools.hxx>
#include "itkNumericSeriesFileNames.h"
#include "itkImageSeriesWriter.h"
#include "itkExceptionObject.h"

typedef itk::Image<float, 2> Image2DType;
typedef itk::Image<float, 3> Image3DType;
typedef itk::ImageSeriesReader<Image3DType> ReaderType;
typedef itk::GDCMImageIO ImageIOType;
typedef itk::GDCMSeriesFileNames NamesGeneratorType;
typedef std::vector<std::string> SeriesIDContainer;
typedef std::vector<std::string> FileNamesContainer;
typedef itk::ThresholdImageFilter<Image3DType> ThresholdImageFilterType;
typedef itk::NumericSeriesFileNames OutputNamesGeneratorType;
typedef itk::ImageSeriesWriter<Image3DType, Image2DType> SeriesWriterType;
#endif

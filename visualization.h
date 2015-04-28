#pragma once

#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"

#include "itkRelabelComponentImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"

#include "itkExtractImageFilter.h"
#include "QuickView.h"
#include "itkCustomColormapFunction.h"
#include "itkScalarToRGBColormapImageFilter.h"
#include "itkRGBPixel.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkThresholdImageFilter.h"

#include <ostream>
#include <random>

#include "verbosityConstant.h"

template< typename TImageType >
static int getSlice(const typename TImageType::Pointer& image, typename itk::Image< typename TImageType::PixelType, 2>::Pointer& slice, const unsigned int sliceNumber) {

  typedef itk::Image< typename TImageType::PixelType, 2> TImage2DType;
  typedef itk::ExtractImageFilter< TImageType, TImage2DType > ExtractFilterType;
  typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();

  typename TImageType::RegionType inputRegion = image->GetLargestPossibleRegion();
  typename TImageType::SizeType size;
  size.Fill(0);
  size[0] = inputRegion.GetSize()[0];
  size[1] = inputRegion.GetSize()[1];


  typename TImageType::IndexType start;
  start.Fill(0);
  start[2] = sliceNumber;

  typename TImageType::RegionType desiredRegion;
  desiredRegion.SetSize(  size  );
  desiredRegion.SetIndex( start );

  extractFilter->SetInput(image);
  extractFilter->SetExtractionRegion( desiredRegion );
  extractFilter->SetDirectionCollapseToIdentity();

  try {
    extractFilter->Update();
  } catch( itk::ExceptionObject & error ) {
    std::cerr << __FILE__ << ":" << __LINE__ << " Error: " << error << std::endl;
    return EXIT_FAILURE;
  }

  slice = extractFilter->GetOutput();

  return EXIT_SUCCESS;
}

template< typename TImage >
int visualize(const typename TImage::Pointer& inputImage) {
  // VISUALIZATION SHIT

  // use this if we don't need more than 255 values on the connected components number of objects
  //typedef ImageType BigImageType;
  const unsigned int Dimension = TImage::ImageDimension;
  typedef itk::Image< int, Dimension > BigImageType;

  typedef itk::RGBPixel<unsigned char> RGBPixelType;
  typedef itk::Image<RGBPixelType, Dimension>  RGBImageType;

  typedef itk::Image<RGBPixelType, 2>  RGBImage2DType;

  typedef itk::Function::CustomColormapFunction<
      typename TImage::PixelType, RGBPixelType > ColormapType;

  typedef float FloatPixelType;
  typedef itk::Image< FloatPixelType, Dimension > FloatImageType;

  typedef itk::Image< FloatPixelType, 2 > FloatImage2DType;

  typedef itk::Image< typename TImage::PixelType, 2> TImage2D;

  /*split the output into several segmentations*/
  //using unsigned char as output image limits us to 255 results
  typedef itk::ConnectedComponentImageFilter < TImage, BigImageType >
      ConnectedComponentImageFilterType;

  typename ConnectedComponentImageFilterType::Pointer connectedFilter = ConnectedComponentImageFilterType::New ();
  connectedFilter->SetInput(inputImage);

  if(VerbosityConstant::verbosity >= VerbosityConstant::HIGH){
    typedef itk::ImageFileWriter< BigImageType > vWriterFilterType;
    typename vWriterFilterType::Pointer vwriter = vWriterFilterType::New();
    vwriter->SetFileName("connected.mha");
    vwriter->SetInput(connectedFilter->GetOutput());
    try {
      vwriter->Update();
    } catch( itk::ExceptionObject & error ) {
      std::cerr << __FILE__ << __LINE__ << "Error: " << error << std::endl;
      return EXIT_FAILURE;
    }
  }


  typedef itk::RelabelComponentImageFilter< BigImageType, BigImageType> RelabelFilterType;
  typename RelabelFilterType::Pointer relabelFilter = RelabelFilterType::New();
  relabelFilter->SetInput(connectedFilter->GetOutput());

  if(VerbosityConstant::verbosity >= VerbosityConstant::HIGH){
    typedef itk::ImageFileWriter< BigImageType > vWriterFilterType;
    typename vWriterFilterType::Pointer vwriter = vWriterFilterType::New();
    vwriter->SetFileName("relabeled.nrrd");
    vwriter->SetInput(relabelFilter->GetOutput());
    try {
      vwriter->Update();
    } catch( itk::ExceptionObject & error ) {
      std::cerr << __FILE__ << __LINE__ << "Error: " << error << std::endl;
      return EXIT_FAILURE;
    }
  }

  typedef itk::RescaleIntensityImageFilter< BigImageType > RescaleFilterType;
  typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput(relabelFilter->GetOutput());
  rescaleFilter->SetOutputMinimum(100);
  rescaleFilter->SetOutputMaximum(255);

  if(VerbosityConstant::verbosity >= VerbosityConstant::HIGH){
    typedef itk::ImageFileWriter< BigImageType > vWriterFilterType;
    typename vWriterFilterType::Pointer vwriter = vWriterFilterType::New();
    vwriter->SetFileName("rescaled.nrrd");
    vwriter->SetInput(rescaleFilter->GetOutput());
    try {
      vwriter->Update();
    } catch( itk::ExceptionObject & error ) {
      std::cerr << __FILE__ << __LINE__ << "Error: " << error << std::endl;
      return EXIT_FAILURE;
    }
  }


  typedef itk::ThresholdImageFilter < BigImageType > ThresholdImageFilterType;

  typename ThresholdImageFilterType::Pointer thresholdFilter
      = ThresholdImageFilterType::New();
  thresholdFilter->SetInput(rescaleFilter->GetOutput());
  thresholdFilter->ThresholdOutside(101, 255);
  thresholdFilter->SetOutsideValue(0);

  typedef itk::ScalarToRGBColormapImageFilter<BigImageType, RGBImageType> ColormapFilterType;
  typename ColormapFilterType::Pointer colormapFilter = ColormapFilterType::New();

  colormapFilter->SetColormap( ColormapFilterType::Hot );
  colormapFilter->SetInput(thresholdFilter->GetOutput());

  try {

    colormapFilter->Update();
    if(VerbosityConstant::verbosity >= VerbosityConstant::LOW)
      std::cout << "Number of objects: " << connectedFilter->GetObjectCount() << std::endl;

  } catch( itk::ExceptionObject & error ) {
    std::cerr << __FILE__ << __LINE__ << "Error: " << error << std::endl;
    return EXIT_FAILURE;
  }

    if(VerbosityConstant::verbosity >= VerbosityConstant::HIGH){
    typedef itk::ImageFileWriter< RGBImageType > vWriterFilterType;
    typename vWriterFilterType::Pointer vwriter = vWriterFilterType::New();
    vwriter->SetFileName("colored.tif");
    vwriter->SetInput(colormapFilter->GetOutput());

    try {
      vwriter->Update();
    } catch( itk::ExceptionObject & error ) {
      std::cerr << __FILE__ << __LINE__ << "Error: " << error << std::endl;
      return EXIT_FAILURE;
    }

  }

  QuickView viewer;

  if(Dimension <= 2) {

    viewer.AddImage(inputImage.GetPointer(), true, "Original");
    viewer.AddRGBImage(colormapFilter->GetOutput(), true, "Relabeled");

  } else {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, inputImage->GetLargestPossibleRegion().GetSize()[2] - 1 );

    const unsigned int sliceNumber = 4;//dis(gen);

    if(VerbosityConstant::verbosity >= VerbosityConstant::MEDIUM)
      std::cout << "Selected slice: " << sliceNumber << std::endl;

    typename TImage2D::Pointer slice1 = TImage2D::New();

    if(getSlice< TImage >(inputImage, slice1, sliceNumber) == EXIT_FAILURE)
      return EXIT_FAILURE;

    viewer.AddImage(slice1.GetPointer(), true, "Original");

    try {

      colormapFilter->Update();

    } catch( itk::ExceptionObject & error ) {
      std::cerr << __FILE__ << ":" << __LINE__ << "Error: " << error << std::endl;
      return EXIT_FAILURE;
    }

    RGBImage2DType::Pointer slice2 = RGBImage2DType::New();

    if(getSlice<RGBImageType>(colormapFilter->GetOutput(), slice2, sliceNumber) == EXIT_FAILURE)
        return EXIT_FAILURE;

    viewer.AddRGBImage(slice2.GetPointer(), true, "Relabeled");

  }

  try {

    viewer.Visualize();

  } catch( itk::ExceptionObject & error ) {
    std::cerr << __FILE__ << ":" << __LINE__ << " Error: " << error << std::endl;
    return EXIT_FAILURE;
  }

}



//////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2015 Pol Monso Purti                                           //
// Contact <pol.monso@gmail.com> for comments & bug reports                   //
//                                                                              //
// This program is free software: you can redistribute it and/or modify         //
// it under the terms of the version 3 of the GNU General Public License        //
// as published by the Free Software Foundation.                                //
//                                                                              //
// This program is distributed in the hope that it will be useful, but          //
// WITHOUT ANY WARRANTY; without even the implied warranty of                   //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU             //
// General Public License for more details.                                     //
//                                                                              //
// You should have received a copy of the GNU General Public License            //
// along with this program. If not, see <http://www.gnu.org/licenses/>.         //
//////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2015 Pol Monso Purti                                           //
// Contact <pol.monso@gmail.com> for comments & bug reports                   //
//                                                                              //
// This program is free software: you can redistribute it and/or modify         //
// it under the terms of the version 3 of the GNU General Public License        //
// as published by the Free Software Foundation.                                //
//                                                                              //
// This program is distributed in the hope that it will be useful, but          //
// WITHOUT ANY WARRANTY; without even the implied warranty of                   //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU             //
// General Public License for more details.                                     //
//                                                                              //
// You should have received a copy of the GNU General Public License            //
// along with this program. If not, see <http://www.gnu.org/licenses/>.         //
//////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "itkImage.h"
#include "itkImageIterator.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"

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

#include "itkRegionOfInterestImageFilter.h"

#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"

#include <ostream>
#include <random>

#include "verbosityConstant.h"

#include "graphCutsAdapter.h"

//we assume that the input is a labelmap
template< typename TImageType >
int GraphCutsAdapter< TImageType >::getBoundingBox(const TImageType* segmentationImage,
                                                   typename TImageType::RegionType& roi,
                                                   typename TImageType::PixelType labelId)
{
  typedef itk::BinaryThresholdImageFilter < TImageType, TImageType > BinaryThresholdImageFilterType;

  typename BinaryThresholdImageFilterType::Pointer binaryThresholdFilter
      = BinaryThresholdImageFilterType::New();
  binaryThresholdFilter->SetInput(segmentationImage);
  binaryThresholdFilter->SetLowerThreshold(labelId);
  binaryThresholdFilter->SetUpperThreshold(labelId);
  binaryThresholdFilter->SetInsideValue(255);
  binaryThresholdFilter->SetOutsideValue(0);

  //TODO do we need to update here or can we wait?
  try {
    binaryThresholdFilter->Update();
  } catch( itk::ExceptionObject & error ) {
    std::cerr << __FILE__ << __LINE__ << "Error: " << error << std::endl;
    return FUCKEDUP;
  }

  if(VerbosityConstant::verbosity >= VerbosityConstant::HIGH){
    typedef itk::ImageFileWriter< TImageType > vWriterFilterType;
    typename vWriterFilterType::Pointer vwriter = vWriterFilterType::New();
    vwriter->SetFileName("thresholded.tif");
    vwriter->SetInput(binaryThresholdFilter->GetOutput());

    try {
      vwriter->Update();
    } catch( itk::ExceptionObject & error ) {
      std::cerr << __FILE__ << __LINE__ << "Error: " << error << std::endl;
      return FUCKEDUP;
    }
  }
  return GraphCutsAdapter< TImageType >::getBoundingBox(binaryThresholdFilter->GetOutput(),
                                                        roi);
  return NAILEDIT;
}

//we assume that there's only one connected component
template< typename TImageType >
int GraphCutsAdapter< TImageType >::getBoundingBox(const TImageType* binarySegmentationImage,
                                                   typename TImageType::RegionType& roi)
{

  typedef itk::BinaryImageToShapeLabelMapFilter< TImageType > BinaryImageToShapeLabelMapFilterType;
  typename BinaryImageToShapeLabelMapFilterType::Pointer binaryImageToShapeLabelMapFilter = BinaryImageToShapeLabelMapFilterType::New();
  binaryImageToShapeLabelMapFilter->SetInput(binarySegmentationImage);

  try {

    binaryImageToShapeLabelMapFilter->Update();

  } catch( itk::ExceptionObject & error ) {
    std::cerr << __FILE__ << ":" << __LINE__ << "Error: " << error << std::endl;
    return FUCKEDUP;
  }

  const unsigned int numObjects = binaryImageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects();

  //we assume that there's only one connected component
  assert(numObjects == 1  && "getBoundingBox assumes you filtered out all other labels");

  if(VerbosityConstant::verbosity >= VerbosityConstant::MEDIUM)
    std::cout << "There are " << numObjects << " objects." << std::endl;

  if(numObjects == 0) {
    std::cerr << "0 objects present. Stop." << std::endl;
    return FUCKEDUP;
  }

  typedef typename BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType LabelObjectType;

  LabelObjectType* labelObject = binaryImageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(0);

  if(VerbosityConstant::verbosity >= VerbosityConstant::MEDIUM) {
    std::cout << "LabelObject description:"        << std::endl;
    std::cout << "Object "        << std::endl;
    std::cout << " Bounding box " << labelObject->GetBoundingBox() << std::endl;
    std::cout << " Centroid "     << labelObject->GetCentroid() << std::endl;
    std::cout << " Volume "       << labelObject->GetPhysicalSize() << std::endl;
  }

  roi = labelObject->GetBoundingBox();

  return NAILEDIT;
}


//we assume binary image
template< typename TImageType >
int GraphCutsAdapter< TImageType >::seedSinksImage2Coordinates(const TImageType* seedSinkROI,
                                                               std::vector< typename TImageType::IndexType >& seeds,
                                                               std::vector< typename TImageType::IndexType >& sinks) {

  //extract cube or assume extracted?
  if(VerbosityConstant::verbosity >= VerbosityConstant::HIGH){
    typename WriterFilterType::Pointer writer = WriterFilterType::New();
    writer->SetFileName("seedSinks.tif");
    writer->SetInput(seedSinkROI);

    try {
      writer->Update();
    } catch( itk::ExceptionObject & error ) {
      std::cerr << __FILE__ << __LINE__ << "Error: " << error << std::endl;
      return FUCKEDUP;
    }
  }

  //retrieve a list of coordinates from seedSink
  itk::ImageRegionConstIteratorWithIndex< TImageType > it(seedSinkROI, seedSinkROI->GetLargestPossibleRegion());
  for(it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    if(it.Get() == SEED_LABEL){
      seeds.push_back( it.GetIndex() );
    } else if(it.Get() == SINK_LABEL){
      sinks.push_back( it.GetIndex() );
    }
  }

  if(sinks.size() == 0 || seeds.size() == 0){
    if(VerbosityConstant::verbosity >= VerbosityConstant::LOW)
      std::cerr << "No seeds or sinks found." << std::endl;
    return FUCKEDUP;
  }
  return NAILEDIT;

}

template< typename TImageType >
int GraphCutsAdapter< TImageType >::bilabelImage2LabelObjects(const TImageType* segmentationImage,
                                                              typename ShapeLabelObjectType::Pointer& labelObject1,
                                                              typename ShapeLabelObjectType::Pointer& labelObject2) {

  typedef itk::LabelImageToShapeLabelMapFilter< TImageType > ImageToShapeLabelMapFilterType;
  typename ImageToShapeLabelMapFilterType::Pointer imageToShapeLabelMapFilter = ImageToShapeLabelMapFilterType::New();
  imageToShapeLabelMapFilter->SetInput(segmentationImage);

  try {

    imageToShapeLabelMapFilter->Update();

  } catch( itk::ExceptionObject & error ) {
    std::cerr << __FILE__ << ":" << __LINE__ << "Error: " << error << std::endl;
    return FUCKEDUP;
  }

  const unsigned int numObjects = imageToShapeLabelMapFilter->GetOutput()->GetNumberOfLabelObjects();

  //we assume that there's only one connected component
  assert(numObjects == 2  && " we assume you only have 2 labels, but there weren't.");

  if(numObjects != 2) {
    std::cerr << numObjects << " !=2 objects present. Stop." << std::endl;
    return JUSTONEOBJECT;
  }

  typedef typename ImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType LabelObjectType;

  labelObject1 = imageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(0);
  labelObject2 = imageToShapeLabelMapFilter->GetOutput()->GetNthLabelObject(1);

  if(VerbosityConstant::verbosity >= VerbosityConstant::MEDIUM) {
    std::cout << "BilabelImage to Label Objects:"        << std::endl;
    std::cout << "Object 1"        << std::endl;
    std::cout << " Bounding box " << labelObject1->GetBoundingBox() << std::endl;
    std::cout << " Centroid "     << labelObject1->GetCentroid() << std::endl;
    std::cout << " Volume "       << labelObject1->GetPhysicalSize() << std::endl;
    std::cout << "Object 2"        << std::endl;
    std::cout << " Bounding box " << labelObject2->GetBoundingBox() << std::endl;
    std::cout << " Centroid "     << labelObject2->GetCentroid() << std::endl;
    std::cout << " Volume "       << labelObject2->GetPhysicalSize() << std::endl;
  }

  return NAILEDIT;

}

//TODO avoid copy of std::vectors
template< typename TImageType >
int GraphCutsAdapter< TImageType >::dummygraphcuts(const TImageType* segmentationImage,
                                                   const std::vector< typename GradientImageType::Pointer >& scoreImages,
                                                   std::vector< typename TImageType::IndexType > remappedseeds,
                                                   std::vector< typename TImageType::IndexType > remappedsinks,
                                                   typename TImageType::Pointer& splittedSegmentationImage) {

  assert(segmentationImage->GetLargestPossibleRegion().GetSize()[0] > 0);

  //TODO maybe it's ok to delegate the allocation to the client
  splittedSegmentationImage->SetRegions(segmentationImage->GetLargestPossibleRegion());
  splittedSegmentationImage->Allocate();

  //Just checking the gradient
  if(VerbosityConstant::verbosity >= VerbosityConstant::HIGH){
    typedef itk::ImageFileWriter< GradientImageType > GradientWriterFilterType;
    typename GradientWriterFilterType::Pointer gwriter = GradientWriterFilterType::New();

    for(unsigned int i = 0; i < scoreImages.size() ; i++) {
      gwriter->SetInput(scoreImages[i]);
      gwriter->SetFileName("gradient" + std::to_string(i) + ".mha");
      try {
        gwriter->Update();
      } catch( itk::ExceptionObject & error ) {
        std::cerr << __FILE__ << ":" << __LINE__ << " Error: " << error << std::endl;
        return FUCKEDUP;
      }
    }
  }

  //we don't have the graphcuts yet, so let's say the result ofthe graphcut is just the sinkSeedImage

  if(VerbosityConstant::verbosity >= VerbosityConstant::HIGH)
    std::cout << "Seeds " << remappedseeds.size() << " Sinks " << remappedsinks.size() << std::endl;

  assert(remappedseeds.size() > 0 && remappedsinks.size() > 0);

  for(int i = 0; i<remappedseeds.size(); i++)
    splittedSegmentationImage->SetPixel(remappedseeds[i], 128);

  for(int i = 0; i<remappedsinks.size(); i++)
    splittedSegmentationImage->SetPixel(remappedsinks[i], 255);

  return NAILEDIT;
}

//template< typename TImageType >
//int GraphCutsAdapter< TImageType >::mergeRegions(std::vector<typename TImageType::RegionType *> regions,
//                                                 typename TImageType::RegionType& mergedRegion)
//{



//}

// Assumes labelMapImage is already allocated (otherwise GetLargestPossibleRegion returns size = 0)
template< typename TImageType >
int GraphCutsAdapter< TImageType >::labelObjects2Image(ShapeLabelObjectType* labelObject1,
                                                       ShapeLabelObjectType* labelObject2,
                                                       typename TImageType::Pointer& labelMapImage){


  //TODO maybe there's a better way to check allocation/size

  typename TImageType::RegionType region = labelMapImage->GetLargestPossibleRegion();
  typename TImageType::SpacingType spacing = labelMapImage->GetSpacing();

  if(region.GetSize()[0] == 0) {
    std::cerr << "Label image not allocated. Using labelobjects boundingbox." << std::endl;
    mergeRegions(std::vector<typename TImageType::RegionType>{labelObject1->GetBoundingBox(), labelObject2->GetBoundingBox()}, region);
  }

  auto segLabelMap = LabelMapType::New();
  segLabelMap->SetSpacing(spacing);
  segLabelMap->SetRegions(region);
  segLabelMap->Allocate();

  labelObject1->SetLabel(itk::NumericTraits< typename TImageType::PixelType >::max());
  labelObject2->SetLabel(itk::NumericTraits< typename TImageType::PixelType >::max()/2);

  segLabelMap->AddLabelObject(labelObject1);
  segLabelMap->AddLabelObject(labelObject2);

  auto label2volume = Label2VolumeFilter::New();
  label2volume->SetInput(segLabelMap);

  try {
    label2volume->Update();
  } catch( itk::ExceptionObject & error ) {
    std::cerr << __FILE__ << ":" << __LINE__ << " Error: " << error << std::endl;
    return FUCKEDUP;
  }
  labelMapImage = label2volume->GetOutput();

  return NAILEDIT;
}

template< typename TImageType >
int GraphCutsAdapter< TImageType >::process(const TImageType* image,
                                            const TImageType* segmentationImage,
                                            const std::vector< typename TImageType::IndexType >& seeds,
                                            const std::vector< typename TImageType::IndexType >& sinks,
                                            typename ShapeLabelObjectType::Pointer& labelObject1,
                                            typename ShapeLabelObjectType::Pointer& labelObject2,
                                            float shapeWeight) {

  //if we allow Real Time seeds might yet not be both there
  if(seeds.size() == 0 || sinks.size() == 0)
    return YOUARENOTREADY;

  //extract cube
  typename TImageType::RegionType roi;
  getBoundingBox(segmentationImage, roi);

  typedef itk::RegionOfInterestImageFilter< TImageType, TImageType > ROIFilterType;
  typename ROIFilterType::Pointer imageROIextractor = ROIFilterType::New();
  imageROIextractor->SetRegionOfInterest(roi);
  imageROIextractor->SetInput(image);

  //features
  typedef itk::RecursiveGaussianImageFilter< TImageType, GradientImageType > GradientFilterType;
  typename GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
  gradientFilter->SetInput(imageROIextractor->GetOutput());
  float sigma = 3.5; //TODO use 3.5, pass by parameter or what?
  gradientFilter->SetSigma( sigma );
  gradientFilter->SetFirstOrder();
  std::vector< typename GradientImageType::Pointer > gradients(TImageType::ImageDimension);

  typedef itk::RegionOfInterestImageFilter< TImageType, TImageType > ROIFilterType;
  typename ROIFilterType::Pointer segmentationROIextractor = ROIFilterType::New();
  segmentationROIextractor->SetRegionOfInterest(roi);
  segmentationROIextractor->SetInput(segmentationImage);

  try {
    segmentationROIextractor->Update();

    // unnecessary if we have pipeline
    // imageROIextractor->Update();
    //TODO load from disk if available
    for(unsigned int dim = 0; dim < TImageType::ImageDimension; dim++) {
      gradientFilter->SetDirection(dim);
      gradientFilter->Update();
      gradients[dim] = gradientFilter->GetOutput();
      gradients[dim]->DisconnectPipeline();
    }
  } catch( itk::ExceptionObject & error ) {
    std::cerr << __FILE__ << ":" << __LINE__ << " Error: " << error << std::endl;
    return FUCKEDUP;
  }

  typename TImageType::Pointer segmentationROI = segmentationROIextractor->GetOutput();

  //remap coordinates
  std::vector< typename TImageType::IndexType > remappedseeds;
  std::vector< typename TImageType::IndexType > remappedsinks;
  remappedseeds.reserve(seeds.size());
  remappedsinks.reserve(sinks.size());

  //FIXME find a way to transform offsets to indeces or a better way
  // to remap the coordinates, this is ridiculous
  // see http://public.kitware.com/pipermail/insight-developers/2006-December/008806.html
  const typename TImageType::IndexType roiBeginIndex = roi.GetIndex();
  for(unsigned int i=0; i < seeds.size(); i++) {
    typename TImageType::OffsetType offset = seeds[i] - roiBeginIndex;
    typename TImageType::IndexType remappedidx;
    for(int dim = 0; dim < TImageType::ImageDimension; dim++)
      remappedidx[dim] = offset[dim];
    remappedseeds.push_back(remappedidx);
  }
  for(unsigned int i=0; i < sinks.size(); i++) {
    typename TImageType::OffsetType offset = sinks[i] - roiBeginIndex;
    typename TImageType::IndexType remappedidx;
    for(int dim = 0; dim < TImageType::ImageDimension; dim++)
      remappedidx[dim] = offset[dim];
    remappedsinks.push_back(remappedidx);
  }

  //normalisation
#define TESTSPACING 0
#if TESTSPACING
  typename TImageType::SpacingType spacing = image->GetSpacing();
  spacing[2] = 5;
#else
  const typename TImageType::SpacingType spacing = image->GetSpacing();
#endif
  //TODO add the possibility to pass spacing as a parameter (as e.g. tif doesn't store it)
  const float zAnisotropyFactor = 2*spacing[2]/(spacing[0] + spacing[1]);

  assert(TImageType::ImageDimension > 0);
  const unsigned int dims = TImageType::ImageDimension;
  for(unsigned int dim = 0; dim < dims; dim++)
    assert(gradients[dim]->GetRequestedRegion().GetSize() == roi.GetSize());

  //compute 1/p_z*(1/(1+|g|) + shapeWeight)
  //we can also use addition filter and inverse filter (si jamais)
  //x and y dimensions
  for(unsigned int dim = 0; dim < dims-1; dim++) {
    itk::ImageRegionIterator< GradientImageType > git( gradients[dim], gradients[dim]->GetRequestedRegion() );
    while( !git.IsAtEnd()) {
        git.Set( 1/(1+std::abs(git.Get())) + shapeWeight );
        ++git;
    }
  }

  //z dimension
  {
    itk::ImageRegionIterator< GradientImageType > git( gradients[dims-1], gradients[dims-1]->GetRequestedRegion() );
    while( !git.IsAtEnd()) {
      git.Set((1/zAnisotropyFactor)*(1/(1+std::abs(git.Get())) + shapeWeight ));
      ++git;
    }
  }

  if(VerbosityConstant::verbosity >= VerbosityConstant::HIGH){

    std::cout << "Gradients computed. Saving." << std::endl;

    for(unsigned int dim = 0; dim < dims; dim++) {

      typename GradientWriterFilterType::Pointer writer = GradientWriterFilterType::New();
      writer->SetInput(gradients[dim]);
      writer->SetFileName("gradient" + std::to_string(dim) + ".tif");
      try {
        writer->Update();
      } catch( itk::ExceptionObject & error ) {
        std::cerr << __FILE__ << ":" << __LINE__ << " Error: " << error << std::endl;
        return FUCKEDUP;
      }
    }
  }

  //Do graphcuts

  typename TImageType::Pointer cutSegmentationImage = TImageType::New();

  //temporary patch
  int result;
  {
    result = dummygraphcuts(segmentationROI, gradients, remappedseeds, remappedsinks, cutSegmentationImage);
  }
  //
  if(result == FUCKEDUP)
    return result;

  if(VerbosityConstant::verbosity >= VerbosityConstant::HIGH){
    typename WriterFilterType::Pointer writer = WriterFilterType::New();
    writer->SetInput(cutSegmentationImage);
    writer->SetFileName("graphCutsOutput.tif");
    try {
      writer->Update();
    } catch( itk::ExceptionObject & error ) {
      std::cerr << __FILE__ << ":" << __LINE__ << " Error: " << error << std::endl;
      return FUCKEDUP;
    }
  }

  //Transform graphcut output to labelObjects
  result = bilabelImage2LabelObjects(cutSegmentationImage, labelObject1, labelObject2);

  if(result == FUCKEDUP || result == JUSTONEOBJECT)
    return result;

  if(VerbosityConstant::verbosity >= VerbosityConstant::HIGH){

    typename TImageType::Pointer labelMapImage = TImageType::New();

    labelObjects2Image(labelObject1, labelObject2, labelMapImage);

    typename WriterFilterType::Pointer writer = WriterFilterType::New();
    writer->SetFileName("ROILabelMap.tif");
    writer->SetInput(labelMapImage);
    try {
      writer->Update();
    } catch( itk::ExceptionObject & error ) {
      std::cerr << __FILE__ << ":" << __LINE__ << " Error: " << error << std::endl;
      return FUCKEDUP;
    }
  }

  //reset the indexes to its original Region
  typename TImageType::OffsetType offset;
  typename TImageType::IndexType index = roi.GetIndex();
  for(unsigned int i = 0; i < TImageType::ImageDimension; i++)
    offset[i] = index[i];

#warning Shift doesn't work correctly given that \
  it does not change the bounding box. \
  Setting the boundingbox does not work either, \
  given that the indexes of the pixels remain unchanged. \
  Therefore, we do both until the bug is solved.

#define USESHIFT
#ifdef USESHIFT
  labelObject1->Shift(offset);
  labelObject2->Shift(offset);

  //we need to do this in order to have updated attributes.
  //a bit of an overkill, so maybe consider returning a labelmap instead of two objects
//  ShapeLabelMapFilter::Pointer shape = ShapeFilterType::New();
//  shape->SetInput(labelMap);
//  shape->Update();
//#else
  typename TImageType::RegionType region1 = labelObject1->GetBoundingBox();
  region1.GetModifiableIndex() += offset;
  labelObject1->SetBoundingBox(region1);
  typename TImageType::RegionType region2 = labelObject2->GetBoundingBox();
  region2.GetModifiableIndex() += offset;
  labelObject2->SetBoundingBox(region2);
#endif
  return NAILEDIT;

}

//explicit instantiation will fail if the caller uses a different instantiation
template class GraphCutsAdapter<itk::Image<unsigned char, 3> >; // explicit template instantiation

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
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
//#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"

#include "itkRegionOfInterestImageFilter.h"

#define YOUARENOTREADY 2
#define JUSTONEOBJECT 1
#define NAILEDIT 0
#define FUCKEDUP -1

#define SEED_LABEL 255
#define SINK_LABEL 128

template< typename TImageType >
class GraphCutsAdapter {
  public:
    typedef itk::Image<float, TImageType::ImageDimension > GradientImageType;
    typedef itk::ImageFileWriter< GradientImageType > GradientWriterFilterType;
    typedef itk::CovariantVector< float, TImageType::ImageDimension > CoVectorType;
    typedef itk::Image< CoVectorType, TImageType::ImageDimension > ChanneledGradientImageType;

    typedef itk::ImageFileWriter< TImageType > WriterFilterType;
    typedef itk::ImageFileReader< TImageType > ReaderFilterType;

    typedef itk::ShapeLabelObject< typename TImageType::PixelType, TImageType::ImageDimension > ShapeLabelObjectType;
    typedef itk::LabelMap< ShapeLabelObjectType > LabelMapType;
    typedef itk::LabelImageToShapeLabelMapFilter< TImageType, LabelMapType > Image2LabelFilter;
    typedef itk::LabelMapToLabelImageFilter< LabelMapType, TImageType > Label2VolumeFilter;

    //we assume that the input is a labelmap
    static int getBoundingBox(const TImageType* segmentationImage,
                              typename TImageType::RegionType& roi,
                              typename TImageType::PixelType labelId);

    //we assume that there's only one connected component
    static int getBoundingBox(const TImageType* binarySegmentationImage,
                              typename TImageType::RegionType& roi);


    //we assume binary image
    static int seedSinksImage2Coordinates(const TImageType* seedSinkROI,
                                          std::vector< typename TImageType::IndexType >& seeds,
                                          std::vector< typename TImageType::IndexType >& sinks);

    static int bilabelImage2LabelObjects(const TImageType* segmentationImage,
                                         typename ShapeLabelObjectType::Pointer& labelObject1,
                                         typename ShapeLabelObjectType::Pointer& labelObject2);

    //TODO avoid copy of std::vectors
    /**
     * @brief dummy graphcuts method
     *
     * Copies the seeds and sinks with different labels to the output.
     * A dummy method that replaces the graphcuts which would return the
     * splitted image
     *
     * @param segmentationImage         image with the image data (usually only the ROI)
     * @param scoreImages               image with the gradient of segmentationImage, components x,y and z.
     * @param seeds                     coordinates of the graphcuts' seeds
     * @param sinks                     coordinates of the graphcuts' sinks
     * @param splittedSegmentationImage output image with the splitted segmentation, that is,
     * two segmentations with two different labels or two connected components separated by
     * a background value pixel (not fully connected)
     * @return -1 on failure, 0 on success
     */
    static int dummygraphcuts(const TImageType* segmentationImage,
                              const std::vector< typename GradientImageType::Pointer >& scoreImages,
                              const std::vector< typename TImageType::IndexType >& sources,
                              const std::vector< typename TImageType::IndexType >& sinks,
                              typename TImageType::Pointer& splittedSegmentationImage);
    static int graphcuts(const TImageType* segmentationImage,
                              const std::vector< typename GradientImageType::Pointer >& scoreImages,
                              const std::vector< typename TImageType::IndexType >& sources,
                              const std::vector< typename TImageType::IndexType >& sinks,
                              typename TImageType::Pointer& splittedSegmentationImage);

    static int mergeRegions(const std::vector< typename TImageType::RegionType >& regions,
                            typename TImageType::RegionType& mergedRegion) {

      typename TImageType::IndexType index    = regions[0].GetIndex();
      typename TImageType::SizeType  size     = regions[0].GetSize();
      typename TImageType::IndexType indexEnd = regions[0].GetIndex() + size;

      for(auto region : regions)
      {
        typename TImageType::IndexType indexAux = region.GetIndex();
        typename TImageType::IndexType indexEndAux = region.GetIndex() + region.GetSize();
        for(unsigned int dim = 0; dim < TImageType::ImageDimension; dim++)
        {
          index[dim]    = std::min(indexAux[dim], index[dim]);
          indexEnd[dim] = std::max(indexEndAux[dim], indexEnd[dim]);
        }
      }

      for(unsigned int dim = 0; dim < TImageType::ImageDimension; dim++)
      {
        size[dim] = indexEnd[dim] - index[dim];
      }

      mergedRegion.SetIndex(index);
      mergedRegion.SetSize(size);

    }

    static int labelObjects2Image(ShapeLabelObjectType* labelObject1,
                                  ShapeLabelObjectType* labelObject2,
                                  typename TImageType::Pointer& labelMapImage);
    /**
     * @param image             image with the original image data
     * @param segmentationImage image with the selected segmentation to split
     * @param seeds             vector of coordinates for the seeds
     * @param sinks             vector of coordinates for the sinks
     * @param labelObject1      output label object
     * @param labelObject1      output label object
     * @param shapeWeight       value added to the weight to favor shape uniformity:
     *                          1/p_z*(1/(1+|g|) + shapeWeight)
     */
    static int process(const TImageType* image,
                       const TImageType* segmentationImage,
                       const std::vector< typename TImageType::IndexType >& seeds,
                       const std::vector< typename TImageType::IndexType >& sinks,
                       typename ShapeLabelObjectType::Pointer& labelObject1,
                       typename ShapeLabelObjectType::Pointer& labelObject2,
                       float shapeWeight = 0);
  };

#ifndef ITK_MANUAL_INSTANTIATION
#include "graphCutsAdapter.hxx"
#endif

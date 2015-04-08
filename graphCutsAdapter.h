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

#define NAILEDIT 0
#define FUCKEDUP -1

#define SEED_LABEL 255
#define SINK_LABEL 128

template< typename TImageType >
class GraphCutsAdapter {
  public:
    typedef itk::Image<float, TImageType::ImageDimension > GradientImageType;

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
    static int dummygraphcuts(const TImageType* segmentationImage,
                              const GradientImageType* gradientImage,
                              const std::vector< typename TImageType::IndexType >& seeds,
                              const std::vector< typename TImageType::IndexType >& sinks,
                              typename TImageType::Pointer& cutSegmentationImage);

    static int mergeRegions(std::vector< typename TImageType::RegionType >& regions,
                            typename TImageType::RegionType& region){

      typename TImageType::IndexType index;
      typename TImageType::SizeType size;

      for(unsigned int dim = 0; dim < TImageType::ImageDimension; dim++)
      {
        for(auto region : regions)
          {
          //[dim] = std::max(region.GetIndex()[dim], )

          }
      }
    }

    static int labelObjects2Image(ShapeLabelObjectType* labelObject1,
                                  ShapeLabelObjectType* labelObject2,
                                  typename TImageType::Pointer& labelMapImage);
    /**
       * @param segmentationImage with the scores
       * @param seedSinksImage with the two labels for seeds and sinks
       * @param bilabelImage label image with the two splitted conected components
       */
    static int process(const TImageType* segmentationImage,
                       std::vector< typename TImageType::IndexType > seeds,
                       std::vector< typename TImageType::IndexType > sinks,
                       typename ShapeLabelObjectType::Pointer& labelObject1,
                       typename ShapeLabelObjectType::Pointer& labelObject2);
  };

#ifndef ITK_MANUAL_INSTANTIATION
#include "graphCutsAdapter.hxx"
#endif

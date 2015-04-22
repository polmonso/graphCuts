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


#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"

#include "graphCutsAdapter.h"

#include "verbosityConstant.h"

#include "tclap/CmdLine.h"

#include "Visualization.h"


int main( int argc, char* argv[] )
{

  const unsigned int Dimension = 3;

  typedef unsigned char                          PixelType;
  typedef itk::Image< PixelType, Dimension >     itkVolumeType;
  typedef itk::ImageFileReader<itkVolumeType>    ReaderFilterType;
  typedef itk::ImageFileWriter<itkVolumeType>    WriterFilterType;

  typedef GraphCutsAdapter<itkVolumeType> GCA;
  typedef GCA::ShapeLabelObjectType ShapeLabelObjectType;
  typedef GCA::LabelMapType         LabelMapType;
  typedef GCA::Image2LabelFilter    Image2LabelFilter;
  typedef GCA::Label2VolumeFilter   Label2VolumeFilter;

  //commmand line options
  std::string imageFilename;
  std::string segmentationFilename;
  std::string seedSinkFilename;
  float zspacing;
  int labelid;

  try
  {
    TCLAP::CmdLine cmd("GraphCuts standalone", ' ', "0.1");
    TCLAP::ValueArg<std::string> imageFilenameArg("i","inputVolume","Input Image Volume whose segmentation refers to",true,"../data/image.tif","string");
    TCLAP::ValueArg<std::string> segmentationFilenameArg("s","inputScoreVolume","Input Score Volume to split",true,"../data/score.tif","string");
    TCLAP::ValueArg<std::string> seedSinkFilenameArg("g","seedsAndSinks","Bilabel Volume with lowest label for seeds and highest for sinks",true,"","string");
    TCLAP::ValueArg<float> zspacingArg("z","zspacing","normalized spacing",false,1.0, "float");
    TCLAP::ValueArg<int> selectedSegmentationLabelArg("l","labelid","label id to split",false,1, "int");
    TCLAP::SwitchArg verbosityArg("v","verbose","Verbosity level on", false);

    cmd.add( imageFilenameArg );
    cmd.add( segmentationFilenameArg );
    cmd.add( seedSinkFilenameArg);
    cmd.add( zspacingArg );
    cmd.add( verbosityArg );
    cmd.add( selectedSegmentationLabelArg );
    cmd.parse( argc, argv );

    imageFilename        = imageFilenameArg.getValue();
    segmentationFilename = segmentationFilenameArg.getValue();
    seedSinkFilename     = seedSinkFilenameArg.getValue();
    zspacing             = zspacingArg.getValue();
    labelid              = selectedSegmentationLabelArg.getValue();

    //Constant singleton init
    if(verbosityArg.getValue())
      VerbosityConstant::verbosity = VerbosityConstant::HIGH;
    else
      VerbosityConstant::verbosity = VerbosityConstant::LOW;

    std::ifstream f(imageFilename.c_str());
    if (!f.good()) {
      std::cerr << "file " << imageFilename << " does not exist" << std::endl;
      return EXIT_FAILURE;
    }
    f.close();

    std::ifstream g(segmentationFilename.c_str());
    if (!g.good()) {
      std::cerr << "file " << segmentationFilename << " does not exist" << std::endl;
      return EXIT_FAILURE;
    }
    g.close();

    std::ifstream h(seedSinkFilename.c_str());
    if (!h.good()) {
      std::cerr << "file " << seedSinkFilename << " does not exist" << std::endl;
      return EXIT_FAILURE;
    }
    h.close();

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    return EXIT_FAILURE;
  }
  //
  if(VerbosityConstant::verbosity >= VerbosityConstant::MEDIUM)
    std::cout << "Filename: " << segmentationFilename << std::endl;

  typename ReaderFilterType::Pointer imageReader = ReaderFilterType::New();
  imageReader->SetFileName(imageFilename);

  typename ReaderFilterType::Pointer segmentationReader = ReaderFilterType::New();
  segmentationReader->SetFileName(segmentationFilename);

  typename ReaderFilterType::Pointer seedSinkReader = ReaderFilterType::New();
  seedSinkReader->SetFileName(seedSinkFilename);

  try {
    imageReader->Update();
    segmentationReader->Update();
    seedSinkReader->Update();
  } catch(itk::ExceptionObject &exp) {
    std::cerr << "exception: " << exp.what() << ". Caught at " << __FILE__ << ":" << __LINE__ << std::endl;
    return EXIT_FAILURE;
  }

  typename itkVolumeType::Pointer image = imageReader->GetOutput();
  typename itkVolumeType::Pointer segmentationImage = segmentationReader->GetOutput();
  typename itkVolumeType::Pointer seedSinkImage = seedSinkReader->GetOutput();

  if(zspacing != 1.0) {
    //Note that if we read a mha file the zspacing will be 1.0 by default (no parameter inputed)
    //so the image spacings will be already correct, we will skip this step
    typename itkVolumeType::SpacingType spacing;
    spacing[0] = 1.0;
    spacing[1] = 1.0;
    spacing[2] = zspacing;
    std::cout << "Setting spacing to " << spacing << std::endl;

    image->SetSpacing(spacing);
    segmentationImage->SetSpacing(spacing);
    seedSinkImage->SetSpacing(spacing);
  }

//#define testindices
#ifdef testindices
  //test indices mixup
  itk::ImageRegionIterator<itkVolumeType> imageIterator(seedSinkImage, seedSinkImage->GetLargestPossibleRegion());

  itkVolumeType::IndexType originalIndex;
  while(!imageIterator.IsAtEnd())
  {
    itkVolumeType::PixelType pix = imageIterator.Get();
    if(pix == 255) {
      originalIndex = imageIterator.GetIndex();
      break;
    }
    ++imageIterator;
  }
  //done indeces check
#endif

  //todo use pointers
  std::vector< typename itkVolumeType::IndexType > seeds;
  std::vector< typename itkVolumeType::IndexType > sinks;

  if(VerbosityConstant::verbosity >= VerbosityConstant::MEDIUM)
    std::cout << "Seed sinks 2 coordinates" << std::endl;

  int result;
  result = GCA::seedSinksImage2Coordinates(seedSinkImage,
                                           seeds,
                                           sinks);

  if(result == FUCKEDUP)
    return EXIT_FAILURE;

  GCA::ShapeLabelObjectType::Pointer labelObject1 = GCA::ShapeLabelObjectType::New();
  GCA::ShapeLabelObjectType::Pointer labelObject2 = GCA::ShapeLabelObjectType::New();

  if(VerbosityConstant::verbosity >= VerbosityConstant::MEDIUM)
    std::cout << "Process" << std::endl;

  result = GraphCutsAdapter<itkVolumeType>::process(image,
                                                    segmentationImage,
                                                    seeds,
                                                    sinks,
                                                    labelObject1,
                                                    labelObject2);

  if(result == FUCKEDUP || result == YOUARENOTREADY)
    return EXIT_FAILURE;

  //FIXME deal with JUSTONEOBJECT better
  if(result == JUSTONEOBJECT) //we only have one labelObject
    labelObject2 = labelObject1;

  // VISUALIZATION
  if(VerbosityConstant::verbosity >= VerbosityConstant::LOW)
    std::cout << "Visualization" << std::endl;

  if(VerbosityConstant::verbosity >= VerbosityConstant::MEDIUM)
    std::cout << "LabelObjects to Image" << std::endl;

  typename itkVolumeType::Pointer volume = itkVolumeType::New();
  volume->SetRegions(segmentationImage->GetLargestPossibleRegion());
  volume->Allocate();

  GCA::labelObjects2Image(labelObject1,
                          labelObject2,
                          volume);

#ifdef testindices
  itk::ImageRegionIterator<itkVolumeType> imageIterator2(volume, volume->GetLargestPossibleRegion());
  itkVolumeType::IndexType newIndex;
  while(!imageIterator2.IsAtEnd())
  {
    if(imageIterator2.Get() == 255) {
      newIndex = imageIterator2.GetIndex();
      break;
    }
    ++imageIterator2;
  }
  std::cout << "Original Index: " << originalIndex << " New Index: " << newIndex << std::endl;
  std::cout << "Obj 1 Bounding box: " << labelObject1->GetBoundingBox() << std::endl;

  if(VerbosityConstant::verbosity >= VerbosityConstant::HIGH){

    WriterFilterType::Pointer writer = WriterFilterType::New();
    writer->SetFileName("input.tif");
    writer->SetInput(seedSinkImage);

    try {
      writer->Update();
    } catch( itk::ExceptionObject & error ) {
      std::cerr << __FILE__ << __LINE__ << "Error: " << error << std::endl;
      return FUCKEDUP;
    }
    writer->SetFileName("output.tif");
    writer->SetInput(volume);

    try {
      writer->Update();
    } catch( itk::ExceptionObject & error ) {
      std::cerr << __FILE__ << __LINE__ << "Error: " << error << std::endl;
      return FUCKEDUP;
    }
  }
#endif

  // visualize< itkVolumeType >(volume);

  if(result == FUCKEDUP)
    return EXIT_FAILURE;

  if(VerbosityConstant::verbosity >= VerbosityConstant::LOW)
    std::cout << "Job's Done. Have a Nice Day!" << std::endl;

  return EXIT_SUCCESS;
}





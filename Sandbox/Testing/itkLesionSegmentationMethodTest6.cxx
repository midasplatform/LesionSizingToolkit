/*=========================================================================

  Program:   Lesion Sizing Toolkit
  Module:    itkLesionSegmentationMethodTest6.cxx

  Copyright (c) Kitware Inc. 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

// The test runs a shape detection level set from user supplied seed points
// and then runs the shape detection level set with the results from the 
// fast marching to get the final segmentation.

#include "itkLesionSegmentationMethod.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSpatialObject.h"
#include "itkSpatialObjectReader.h"
#include "itkImageMaskSpatialObject.h"
#include "itkLungWallFeatureGenerator.h"
#include "itkSatoVesselnessSigmoidFeatureGenerator.h"
#include "itkSigmoidFeatureGenerator.h"
#include "itkFastMarchingAndShapeDetectionLevelSetSegmentationModule.h"
#include "itkMinimumFeatureAggregator.h"

int main( int argc, char * argv [] )
{

  if( argc < 3 )
    {
    std::cerr << "Applies fast marhching followed by shape detection. Missing Arguments" << std::endl;
    std::cerr << argv[0] << "\n\tlandmarksFile\n\tinputImage\n\toutputImage ";
    std::cerr << "\n\t[RMSErrorForShapeDetection]"
              << "\n\t[IterationsForShapeDetection]"
              << "\n\t[CurvatureScalingForShapeDetection]"
              << "\n\t[PropagationScalingForShapeDetection]";
    std::cerr << "\n\tstopping time for fast marching";
    std::cerr << "\n\tdistance from seeds for fast marching" << std::endl;
    return EXIT_FAILURE;
    }


  const unsigned int Dimension = 3;
  typedef signed short   InputPixelType;

  typedef itk::Image< InputPixelType, Dimension > InputImageType;

  typedef itk::ImageFileReader< InputImageType > InputImageReaderType;
  InputImageReaderType::Pointer inputImageReader = InputImageReaderType::New();

  inputImageReader->SetFileName( argv[2] );

  try 
    {
    inputImageReader->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }


  typedef itk::LesionSegmentationMethod< Dimension >   MethodType;

  MethodType::Pointer  lesionSegmentationMethod = MethodType::New();
  
  typedef itk::ImageMaskSpatialObject< Dimension > ImageMaskSpatialObjectType;

  ImageMaskSpatialObjectType::Pointer regionOfInterest = ImageMaskSpatialObjectType::New();

  lesionSegmentationMethod->SetRegionOfInterest( regionOfInterest );

  typedef itk::SatoVesselnessSigmoidFeatureGenerator< Dimension > VesselnessGeneratorType;
  VesselnessGeneratorType::Pointer vesselnessGenerator = VesselnessGeneratorType::New();

  typedef itk::LungWallFeatureGenerator< Dimension > LungWallGeneratorType;
  LungWallGeneratorType::Pointer lungWallGenerator = LungWallGeneratorType::New();

  typedef itk::SigmoidFeatureGenerator< Dimension >   SigmoidFeatureGeneratorType;
  SigmoidFeatureGeneratorType::Pointer  sigmoidGenerator = SigmoidFeatureGeneratorType::New();
 
  typedef itk::MinimumFeatureAggregator< Dimension >   FeatureAggregatorType;
  FeatureAggregatorType::Pointer featureAggregator = FeatureAggregatorType::New();
  featureAggregator->AddFeatureGenerator( lungWallGenerator );
  featureAggregator->AddFeatureGenerator( vesselnessGenerator );
  featureAggregator->AddFeatureGenerator( sigmoidGenerator );
  lesionSegmentationMethod->AddFeatureGenerator( featureAggregator );

  typedef MethodType::SpatialObjectType    SpatialObjectType;
  typedef itk::ImageSpatialObject< Dimension, InputPixelType  > InputImageSpatialObjectType;
  InputImageSpatialObjectType::Pointer inputObject = InputImageSpatialObjectType::New();

  InputImageType::Pointer inputImage = inputImageReader->GetOutput();

  inputImage->DisconnectPipeline();

  inputObject->SetImage( inputImage );

  lungWallGenerator->SetInput( inputObject );
  vesselnessGenerator->SetInput( inputObject );
  sigmoidGenerator->SetInput( inputObject );
  lungWallGenerator->SetLungThreshold( -400 );
  vesselnessGenerator->SetSigma( 1.0 );
  vesselnessGenerator->SetAlpha1( 0.5 );
  vesselnessGenerator->SetAlpha2( 2.0 );
  sigmoidGenerator->SetAlpha(  1.0  );
  sigmoidGenerator->SetBeta( -200.0 );
 
  typedef itk::FastMarchingAndShapeDetectionLevelSetSegmentationModule< Dimension > SegmentationModuleType;
  SegmentationModuleType::Pointer  segmentationModule = SegmentationModuleType::New();
  segmentationModule->SetMaximumRMSError( argc > 4 ? atof(argv[4]) : 0.0002 );
  segmentationModule->SetMaximumNumberOfIterations( argc > 5 ? atoi(argv[5]) : 300 );
  segmentationModule->SetCurvatureScaling( argc > 6 ? atof(argv[6]) : 1.0 );
  segmentationModule->SetPropagationScaling( argc > 7 ? atof(argv[7]) : 10.0 );
  segmentationModule->SetStoppingValue( argc > 8 ? atof(argv[8]) : 5.0 );
  segmentationModule->SetDistanceFromSeeds( argc > 9 ? atof(argv[9]) : 2.0 );
  lesionSegmentationMethod->SetSegmentationModule( segmentationModule );

  typedef itk::SpatialObjectReader< 3, unsigned short > SpatialObjectReaderType;
  SpatialObjectReaderType::Pointer landmarkPointsReader = SpatialObjectReaderType::New();
  landmarkPointsReader->SetFileName( argv[1] );
  landmarkPointsReader->Update();
  SpatialObjectReaderType::ScenePointer scene = landmarkPointsReader->GetScene();
  if( !scene )
    {
    std::cerr << "No Scene : [FAILED]" << std::endl;
    return EXIT_FAILURE;
    }
  std::cout << "Number of object in the scene:" << scene->GetNumberOfObjects(1) << std::endl;

  typedef SpatialObjectReaderType::SceneType::ObjectListType     ObjectListType;
  ObjectListType * sceneChildren = scene->GetObjects(999999);
  ObjectListType::const_iterator spatialObjectItr = sceneChildren->begin();

  typedef SegmentationModuleType::InputSpatialObjectType  InputSpatialObjectType; 
  typedef itk::LandmarkSpatialObject< Dimension > LandmarkSpatialObjectType;
  const LandmarkSpatialObjectType * landmarkSpatialObject = NULL;
  while( spatialObjectItr != sceneChildren->end() ) 
    {
    std::string objectName = (*spatialObjectItr)->GetTypeName();
    if( objectName == "LandmarkSpatialObject" )
      {
      landmarkSpatialObject = 
        dynamic_cast< const LandmarkSpatialObjectType * >( spatialObjectItr->GetPointer() );
      }
    spatialObjectItr++;
    }
 
  lesionSegmentationMethod->SetInitialSegmentation( landmarkSpatialObject );
  lesionSegmentationMethod->Update();

  
  typedef SegmentationModuleType::SpatialObjectType           SpatialObjectType;
  typedef SegmentationModuleType::OutputSpatialObjectType     OutputSpatialObjectType;
  typedef SegmentationModuleType::OutputImageType             OutputImageType;

  SpatialObjectType::ConstPointer segmentation = segmentationModule->GetOutput();

  OutputSpatialObjectType::ConstPointer outputObject = 
    dynamic_cast< const OutputSpatialObjectType * >( segmentation.GetPointer() );
  OutputImageType::ConstPointer outputImage = outputObject->GetImage();

  typedef itk::ImageFileWriter< OutputImageType >      OutputWriterType;
  OutputWriterType::Pointer writer = OutputWriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( outputImage );
  writer->UseCompressionOn();

  try 
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

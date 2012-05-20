/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMImageIOFactory.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkMetaImageIOFactory.h"
#include "itkImageSeriesReader.h"
#include "itkEventObject.h"
#include "itkImageToVTKImageFilter.h"
#include "itkMetaImageIOFactory.h"

#include "gdcmUIDGenerator.h"

#include "vtkMassProperties.h"
#include "vtkImageData.h"
#include "vtkMarchingCubes.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkSTLWriter.h"
#include "vtkRenderer.h"
#include "vtkCamera.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkActor.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkImageData.h"
#include "vtkProperty.h"
#include "vtkImagePlaneWidget.h"
#include "vtkCellPicker.h"
#include "vtkPolyDataNormals.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkOutlineSource.h"
#include "vtkCommand.h"
#include "vtkWindowToImageFilter.h"
#include "vtkPNGWriter.h"
#include "itkOrientImageFilter.h"

// This needs to come after the other includes to prevent the global definitions
// of PixelType to be shadowed by other declarations.
#include "itkLesionSegmentationImageFilter8.h"
#include "itkLesionSegmentationCommandLineProgressReporter.h"
#include "LesionSegmentationCLI.h"
#include "LesionSegmentationConfig.h"

#ifdef LSTK_USE_AIM
  #include "itkImageToAIMXMLFilter.h"
#endif

#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

// Contants
const std::string  PatientNameTag = "0010|0010";
const std::string  PatientIdTag = "0010|0020";
const std::string  PatientSexTag = "0010|0040";
const std::string  StudyInstanceUIDTag = "0020|000d";
const std::string  SeriesInstanceUIDTag = "0020|000e";
const std::string  SOPClassUIDTag = "0008|0016";
const std::string  SOPInstanceUIDTag = "0008|0018";
const std::string  SeriesRestriction = "0008|0021";
const float        ContourValue = -0.5;
const unsigned int ImageDimension = LesionSegmentationCLI::ImageDimension;


// typdefs
typedef LesionSegmentationCLI::InputImageType        InputImageType;
typedef LesionSegmentationCLI::RealImageType         RealImageType;
typedef itk::ImageFileReader< InputImageType >       InputReaderType;
typedef itk::ImageFileWriter< RealImageType >        OutputWriterType;
typedef InputImageType                               ImageType;
typedef itk::LandmarkSpatialObject< ImageDimension > LandmarkType;
typedef LandmarkType::PointType                      PointType;
typedef LandmarkType::PointListType                  PointListType;
typedef itk::ImageSeriesReader< ImageType >          ReaderType;
typedef itk::Image< float, ImageDimension >          RealImageType;
typedef itk::GDCMImageIO                             ImageIOType;
typedef itk::GDCMSeriesFileNames                     NamesGeneratorType;
typedef std::vector< std::string >                   SeriesIdContainer;
typedef std::vector< std::string >                   FileNamesContainer;
typedef std::vector< std::string >                   UIDStorageVector;
typedef itk::MetaDataObject< std::string >           MetaDataStringType;
typedef const MetaDataStringType*                    MetaDataStringConstPointer;
typedef const itk::MetaDataObjectBase*               MetaDataUncastedPointer;

typedef itk::LesionSegmentationImageFilter8< InputImageType, RealImageType >
                                                     SegmentationFilterType;


// --------------------------------------------------------------------------
class SwitchVisibilityCallback : public vtkCommand
{
public:
  static SwitchVisibilityCallback *New()
    { return new SwitchVisibilityCallback; }

  void SetActor(vtkActor *aActor)
    {
    this->Actor=aActor;
    }
  void SetRenderWindow(vtkRenderWindow *aRenWin)
    {
    this->RenWin=aRenWin;
    }

  virtual void Execute(vtkObject *vtkNotUsed(caller), unsigned long, void*)
    {
    this->Actor->SetVisibility(1-this->Actor->GetVisibility());
    this->RenWin->Render();
    }

protected:
  vtkActor *Actor;
  vtkRenderWindow *RenWin;
};

// --------------------------------------------------------------------------
class ImageAndMetaDataContainer
{
public:
  ImageAndMetaDataContainer(std::string pName,
                            std::string pId,
                            std::string pSex,
                            std::string stUID,
                            std::string seUID,
                            UIDStorageVector sopClassVect,
                            UIDStorageVector sopInstanceVect,
                            InputImageType::Pointer img) :
    patientName(pName),
    patientId(pId),
    patientSex(pSex),
    studyInstanceUID(stUID),
    seriesInstanceUID(seUID),
    sopClassUIDVector(sopClassVect),
    sopInstanceUIDVector(sopInstanceVect)
  {
    image = img;
    gdcm::UIDGenerator gen;
    currentUID = std::string(gen.Generate());
    std::cout << "Creating with UID: " << currentUID << std::endl;
  }
  
  std::string patientName;
  std::string patientId;
  std::string patientSex;
  std::string studyInstanceUID;
  std::string seriesInstanceUID;
  std::string currentUID;

  UIDStorageVector sopClassUIDVector;
  UIDStorageVector sopInstanceUIDVector;

  LesionSegmentationCLI::InputImageType::Pointer image;
};

// --------------------------------------------------------------------------
ImageAndMetaDataContainer* GetImageAndMetaData( std::string dir, bool ignoreDirection )
{
  ReaderType::Pointer reader = ReaderType::New();
  ImageIOType::Pointer dicomIO = ImageIOType::New();

  reader->SetImageIO( dicomIO );

  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

  nameGenerator->SetUseSeriesDetails( true );
  nameGenerator->AddSeriesRestriction("0008|0021" );
  nameGenerator->SetDirectory( dir );

  try
    {
    std::cout << std::endl << "The directory: " << std::endl;
    std::cout << std::endl << dir << std::endl << std::endl;
    std::cout << "Contains the following DICOM Series: ";
    std::cout << std::endl << std::endl;

    const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();

    SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
    SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
    while( seriesItr != seriesEnd )
      {
      std::cout << seriesItr->c_str() << std::endl;
      seriesItr++;
      }

    std::string seriesIdentifier;
    seriesIdentifier = seriesUID.begin()->c_str();

    std::cout << std::endl << std::endl;
    std::cout << "Now reading series: " << std::endl << std::endl;
    std::cout << seriesIdentifier << std::endl;
    std::cout << std::endl << std::endl;

    FileNamesContainer fileNames;

    fileNames = nameGenerator->GetFileNames( seriesIdentifier );
    reader->SetFileNames( fileNames );

    try
      {
      reader->Update();
      }
    catch (itk::ExceptionObject &ex)
      {
      std::cout << ex << std::endl;
      return NULL;
      }

    ReaderType::DictionaryArrayRawPointer dictOfDicts =
      reader->GetMetaDataDictionaryArray();
    UIDStorageVector sopClassUIDVector;
    sopClassUIDVector.reserve(dictOfDicts->size());
    UIDStorageVector sopInstanceUIDVector;
    sopInstanceUIDVector.reserve(dictOfDicts->size());
    ReaderType::DictionaryArrayType::const_iterator itr;
    for( itr = dictOfDicts->begin(); itr != dictOfDicts->end(); ++itr )
      {
      ReaderType::DictionaryRawPointer dict = *itr;

      MetaDataUncastedPointer sopClassUIDBase = dict->Get(SOPClassUIDTag);
      MetaDataUncastedPointer sopInstanceUIDBase = dict->Get(SOPInstanceUIDTag);

      MetaDataStringConstPointer elClassUID;
      MetaDataStringConstPointer elInstanceUID;
      elClassUID = dynamic_cast<MetaDataStringConstPointer>(sopClassUIDBase);
      elInstanceUID = dynamic_cast<MetaDataStringConstPointer>(sopInstanceUIDBase);

      if( elClassUID )
        {
        sopClassUIDVector.push_back(elClassUID->GetMetaDataObjectValue());
        }
      if( elInstanceUID )
        {
        sopInstanceUIDVector.push_back(elInstanceUID->GetMetaDataObjectValue());
        }
      }

    std::string patientName;
    std::string patientId;
    std::string patientSex;
    std::string studyInstanceUID;
    std::string seriesInstanceUID;
    
    dicomIO->GetValueFromTag( PatientNameTag, patientName );
    dicomIO->GetValueFromTag( PatientIdTag, patientName );
    dicomIO->GetValueFromTag( PatientSexTag, patientSex );
    dicomIO->GetValueFromTag( StudyInstanceUIDTag, studyInstanceUID );
    dicomIO->GetValueFromTag( SeriesInstanceUIDTag, seriesInstanceUID );

    std::cout << "Patient Name: " << patientName << std::endl
              << "Patient ID: " << patientId << std::endl
              << "Patient Sex: " << patientSex << std::endl
              << "Study Instance UID: " << studyInstanceUID << std::endl
              << "Series Instance UID: " << seriesInstanceUID << std::endl;

    ImageType::Pointer image = reader->GetOutput();
    ImageType::DirectionType direction;
    direction.SetIdentity();
    image->DisconnectPipeline();

    if (ignoreDirection)
      {
      std::cout << "Ignoring the direction of the DICOM image and using "
                << "identity." << std::endl;
      image->SetDirection(direction);
      }
    ImageAndMetaDataContainer* output;
    output = new ImageAndMetaDataContainer( patientName,
                                            patientId,
                                            patientSex,
                                            studyInstanceUID,
                                            seriesInstanceUID,
                                            sopClassUIDVector,
                                            sopInstanceUIDVector,
                                            image );
    return output;
    }
  catch (itk::ExceptionObject &ex)
    {
    std::cout << ex << std::endl;
    return NULL;
    }

  return NULL;
}

// --------------------------------------------------------------------------
int ViewImageAndSegmentationSurface(
    LesionSegmentationCLI::InputImageType::Pointer image, vtkPolyData *pd,
    double *roi, LesionSegmentationCLI &args )
{

  std::cout << "Setting up visualization..." << std::endl;


  typedef LesionSegmentationCLI::InputImageType InputImageType;

  typedef itk::ImageToVTKImageFilter< LesionSegmentationCLI::InputImageType > RealITKToVTKFilterType;
  RealITKToVTKFilterType::Pointer itk2vtko = RealITKToVTKFilterType::New();
  itk2vtko->SetInput( image);
  itk2vtko->Update();

  // display the results.
  VTK_CREATE( vtkRenderer, renderer );
  VTK_CREATE( vtkRenderWindow, renWin );
  VTK_CREATE( vtkRenderWindowInteractor, iren );

  renderer->GetActiveCamera()->ParallelProjectionOn();

  renWin->SetSize(600, 600);
  renWin->AddRenderer(renderer);
  iren->SetRenderWindow(renWin);

  // use cell picker for interacting with the image orthogonal views.
  //
  VTK_CREATE(vtkCellPicker , picker);
  picker->SetTolerance(0.005);


  //assign default props to the ipw's texture plane actor
  VTK_CREATE(vtkProperty , ipwProp);


  // Create 3 orthogonal view using the ImagePlaneWidget
  //
  vtkSmartPointer< vtkImagePlaneWidget > imagePlaneWidget[3];
  for (unsigned int i = 0; i < 3; i++)
    {
    imagePlaneWidget[i] = vtkSmartPointer< vtkImagePlaneWidget >::New();

    imagePlaneWidget[i]->DisplayTextOn();
#if VTK_MAJOR_VERSION <= 5
    imagePlaneWidget[i]->SetInput(itk2vtko->GetOutput());
#else
    imagePlaneWidget[i]->SetInputData(itk2vtko->GetOutput());
#endif


    imagePlaneWidget[i]->SetPlaneOrientation(i);
    imagePlaneWidget[i]->SetSlicePosition(pd->GetCenter()[i]);
    imagePlaneWidget[i]->SetPicker(picker);
    imagePlaneWidget[i]->RestrictPlaneToVolumeOn();
    double color[3] = {0,0,0};
    color[i] = 1;
    imagePlaneWidget[i]->GetPlaneProperty()->SetColor(color);
    imagePlaneWidget[i]->SetTexturePlaneProperty(ipwProp);
    imagePlaneWidget[i]->SetResliceInterpolateToLinear();
    imagePlaneWidget[i]->SetWindowLevel(1700,-500);
    imagePlaneWidget[i]->SetInteractor( iren );
    imagePlaneWidget[i]->On();
    }

  imagePlaneWidget[0]->SetKeyPressActivationValue('x');
  imagePlaneWidget[1]->SetKeyPressActivationValue('y');
  imagePlaneWidget[2]->SetKeyPressActivationValue('z');
  

  // Set the background to something grayish
  renderer->SetBackground(0.4392, 0.5020, 0.5647);

  VTK_CREATE(vtkPolyDataMapper , polyMapper);
  VTK_CREATE(vtkActor          , polyActor );

  std::cout << "Computing surface normals for shading.." << std::endl;
  VTK_CREATE( vtkPolyDataNormals, normals );
#if VTK_MAJOR_VERSION <=5
  normals->SetInput(pd);
#else
  normals->SetInputData(pd);
#endif
  normals->SetFeatureAngle(60.0);
  normals->Update();

  polyActor->SetMapper( polyMapper );
#if VTK_MAJOR_VERSION <=5
  polyMapper->SetInput( normals->GetOutput() );
#else
  polyMapper->SetInputData( normals->GetOutput() );
#endif

  VTK_CREATE(vtkProperty , property);
  property->SetAmbient(0.1);
  property->SetDiffuse(0.1);
  property->SetSpecular(0.5);
  property->SetColor(0,1,0);
  property->SetLineWidth(2.0);
  property->SetRepresentationToSurface();

  polyActor->SetProperty( property );
  renderer->AddActor( polyActor );

  if (args.GetOptionWasSet("Wireframe"))
    {
    property->SetRepresentationToWireframe();
    }

  // Bring up the render window and begin interaction.
  renderer->ResetCamera();
  renWin->Render();

  VTK_CREATE( vtkInteractorStyleTrackballCamera, style );
  iren->SetInteractorStyle(style);

  SwitchVisibilityCallback *switchVisibility=SwitchVisibilityCallback::New();
  switchVisibility->SetRenderWindow(renWin);
  switchVisibility->SetActor(polyActor);
  iren->AddObserver(vtkCommand::UserEvent,switchVisibility);
  switchVisibility->Delete();

  if (args.GetOptionWasSet("ShowBoundingBox"))
    {
    VTK_CREATE( vtkOutlineSource, outline );
    outline->SetBounds(roi);
    VTK_CREATE( vtkPolyDataMapper, outlineMapper );
#if VTK_MAJOR_VERSION <= 5
    outlineMapper->SetInput(outline->GetOutput());
#else
    outlineMapper->SetInputData(outline->GetOutput());
#endif
    VTK_CREATE( vtkActor, outlineActor );
    outlineActor->SetMapper(outlineMapper);
    outlineActor->GetProperty()->SetColor(0,1,0);
    renderer->AddActor(outlineActor); 
    }

  std::cout << "Bringing up visualization.." << std::endl;

  if (!args.GetValueAsString("Screenshot").empty())
    {

    double camPos[3][3] = { {1,0,0},{0,-1,0},{0,0,-1} };
    double viewUp[3][3] = { {0,0,1},{0,0,1},{0,-1,0} };
    
    for (unsigned int i = 0; i < 3; i++)
      {
      std::ostringstream os;
      os << args.GetValueAsString("Screenshot") 
         << "_" << i << ".png" << std::ends;
      
      renderer->GetActiveCamera()->SetPosition(camPos[i]);
      renderer->GetActiveCamera()->SetFocalPoint(0,0,0);
      renderer->GetActiveCamera()->SetViewUp(viewUp[i]);

      for (unsigned int j = 0; j < 3; j++)
        {
        imagePlaneWidget[j]->Off();
        }
      imagePlaneWidget[i]->On();
      imagePlaneWidget[i]->SetSlicePosition(pd->GetCenter()[i]);
      renderer->ResetCamera();
      renderer->ResetCameraClippingRange();

      // Reset the camera to the full size of the view for the screenshot
      double parallelScale = 0, bounds[6], l2norm = 0;
      itk2vtko->GetOutput()->GetBounds(bounds);
      for (unsigned int k = 0; k < 3; k++)
        {
        if (k!=i)
          {
          l2norm += ((bounds[2*k+1]-bounds[2*k]) * (bounds[2*k+1]-bounds[2*k]));
          parallelScale = std::max(parallelScale, bounds[2*k+1]-bounds[2*k]);
          }
        }
      renderer->GetActiveCamera()->Zoom(sqrt(l2norm)/parallelScale);

      renWin->Render();

      VTK_CREATE( vtkWindowToImageFilter, w2f );
      w2f->SetInput(renWin);

      VTK_CREATE( vtkPNGWriter, screenshotWriter );

      screenshotWriter->SetFileName(os.str().c_str());
      std::cout << "Screenshot saved to " <<os.str().c_str()<< std::endl;
      screenshotWriter->SetInputConnection(w2f->GetOutputPort());
      screenshotWriter->Write();
      }
    }
  else
    {
    iren->Start();
    }

  return EXIT_SUCCESS;
}

#ifdef LSTK_USE_AIM
void WriteAIM( const std::string& fileName, double contourValue,
               SegmentationFilterType::Pointer seg,
               ImageAndMetaDataContainer* data,
               double volume )
{

  typedef itk::ImageToAIMXMLFilter<RealImageType, InputImageType >
    AIMFilterType;
  std::cout << "Writing the output segmented level set in AIM XML format. "
            << fileName << ". "
            << "The segmentation is an isosurface of this image at a "
            << "value of " << contourValue << "." << std::endl;
  AIMFilterType::Pointer aimFilter = AIMFilterType::New();
  aimFilter->SetContourThreshold( contourValue );
  aimFilter->SetInput( seg->GetOutput() );
  aimFilter->SetReference( data->image );
  aimFilter->SetPatientName( data->patientName );
  aimFilter->SetPatientId( data->patientId );
  aimFilter->SetPatientSex( data->patientSex );
  aimFilter->SetStudyInstanceUID( data->studyInstanceUID );
  aimFilter->SetSeriesInstanceUID( data->seriesInstanceUID );
  aimFilter->SetCurrentUID( data->currentUID );
  aimFilter->SetSOPClassUIDs( data->sopClassUIDVector );
  aimFilter->SetSOPInstanceUIDs( data->sopInstanceUIDVector );
  aimFilter->SetVolume( volume );
  aimFilter->Update();
  
  std::ofstream aimFile;
  aimFile.open( fileName.c_str() );
  aimFile << aimFilter->GetOutput();
  aimFile.close();
}
#endif

// --------------------------------------------------------------------------
int main( int argc, char * argv[] )
{
  //register DICOM and META IO factories
  itk::ObjectFactoryBase::RegisterFactory(itk::GDCMImageIOFactory::New());
  itk::ObjectFactoryBase::RegisterFactory(itk::MetaImageIOFactory::New());

  LesionSegmentationCLI args( argc, argv );

  // Read the volume
  InputReaderType::Pointer reader = InputReaderType::New();
  InputImageType::Pointer image;
  ImageAndMetaDataContainer* data = NULL;

  std::cout << "Reading " << args.GetValueAsString("InputImage") << ".." << std::endl;
  if (!args.GetValueAsString("InputDICOMDir").empty())
    {
    std::cout << "Reading from DICOM dir " << args.GetValueAsString("InputDICOMDir") << ".." << std::endl;
    data = GetImageAndMetaData(
      args.GetValueAsString("InputDICOMDir"),
      args.GetValueAsBool("IgnoreDirection"));
    image = data->image;

    if (!image)
      {
      std::cerr << "Failed to read the input image" << std::endl;
      return EXIT_FAILURE;
      }
    }

  if (!args.GetValueAsString("InputImage").empty())
    {
    reader->SetFileName(args.GetValueAsString("InputImage"));
    reader->Update();
    image = reader->GetOutput();
    }

  //To make sure the tumor polydata aligns with the image volume during
  //vtk rendering in ViewImageAndSegmentationSurface(),
  //reorient image so that the direction matrix is an identity matrix.
  typedef itk::ImageFileReader< InputImageType > InputReaderType;
  itk::OrientImageFilter<InputImageType,InputImageType>::Pointer orienter =
  itk::OrientImageFilter<InputImageType,InputImageType>::New();
  orienter->UseImageDirectionOn();
  InputImageType::DirectionType direction;
  direction.SetIdentity();
  orienter->SetDesiredCoordinateDirection (direction);
  orienter->SetInput(image);
  orienter->Update();
  image = orienter->GetOutput();

  // Set the image object on the args
  args.SetImage( image );


  // Compute the ROI region

  double *roi = args.GetROI();
  InputImageType::RegionType region = image->GetLargestPossibleRegion();

  // convert bounds into region indices
  InputImageType::PointType p1, p2;
  InputImageType::IndexType pi1, pi2;
  InputImageType::IndexType startIndex;
  for (unsigned int i = 0; i < ImageDimension; i++)
    {
    p1[i] = roi[2*i];
    p2[i] = roi[2*i+1];
    }

  image->TransformPhysicalPointToIndex(p1, pi1);
  image->TransformPhysicalPointToIndex(p2, pi2);

  InputImageType::SizeType roiSize;
  for (unsigned int i = 0; i < ImageDimension; i++)
    {
    roiSize[i] = fabs(pi2[i] - pi1[i]);
    startIndex[i] = (pi1[i]<pi2[i])?pi1[i]:pi2[i];
    }
  InputImageType::RegionType roiRegion( startIndex, roiSize );
  std::cout << "ROI region is " << roiRegion << std::endl;
  if (!roiRegion.Crop(image->GetBufferedRegion()))
    {
    std::cerr << "ROI region has no overlap with the image region of"
              << image->GetBufferedRegion() << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "ROI region is " << roiRegion <<  " : covers voxels = "
    << roiRegion.GetNumberOfPixels() << " : " <<
    image->GetSpacing()[0]*image->GetSpacing()[1]*image->GetSpacing()[2]*roiRegion.GetNumberOfPixels()
   << " mm^3" << std::endl;


  // Write ROI if requested
  if (args.GetOptionWasSet("OutputROI"))
    {
    typedef itk::RegionOfInterestImageFilter<
      InputImageType, InputImageType > ROIFilterType;

    ROIFilterType::Pointer roiFilter = ROIFilterType::New();
    roiFilter->SetRegionOfInterest( roiRegion );

    typedef itk::ImageFileWriter< InputImageType > ROIWriterType;
    ROIWriterType::Pointer roiWriter = ROIWriterType::New();
    roiWriter->SetFileName( args.GetValueAsString("OutputROI") );
    roiFilter->SetInput( image );
    roiWriter->SetInput( roiFilter->GetOutput() );

    std::cout << "Writing the ROI image as: " <<
      args.GetValueAsString("OutputROI") << std::endl;
    try
      {
      roiWriter->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "ExceptionObject caught !" << err << std::endl;
      return EXIT_FAILURE;
      }
    }


  // Progress reporting
  typedef itk::LesionSegmentationCommandLineProgressReporter ProgressReporterType;
  ProgressReporterType::Pointer progressCommand =
    ProgressReporterType::New();

  // Run the segmentation filter. Clock ticking...

  std::cout << "\n Running the segmentation filter." << std::endl;
  SegmentationFilterType::Pointer seg = SegmentationFilterType::New();
  seg->SetInput(image);
  seg->SetSeeds(args.GetSeeds());
  seg->SetRegionOfInterest(roiRegion);
  seg->AddObserver( itk::ProgressEvent(), progressCommand );
  if (args.GetOptionWasSet("Sigma"))
    {
    seg->SetSigma(args.GetSigmas());
    }
  seg->SetSigmoidBeta(args.GetValueAsBool("PartSolid") ? -500 : -200 );
  seg->Update();

  if (!args.GetValueAsString("OutputImage").empty())
    {
    std::cout << "Writing the output segmented level set."
      << args.GetValueAsString("OutputImage") <<
      ". The segmentation is an isosurface of this image at a value of -0.5"
      << std::endl;
    OutputWriterType::Pointer writer = OutputWriterType::New();
    writer->SetFileName(args.GetValueAsString("OutputImage"));
    writer->SetInput(seg->GetOutput());
    writer->Update();
    }

  // Compute volume
  typedef itk::ImageToVTKImageFilter< RealImageType > RealITKToVTKFilterType;
  RealITKToVTKFilterType::Pointer itk2vtko = RealITKToVTKFilterType::New();
  itk2vtko->SetInput( seg->GetOutput() );
  itk2vtko->Update();

  std::cout << "Generating an isosurface of the zero-level set (iso-value of -0.5)" << std::endl;
  vtkSmartPointer< vtkMarchingCubes > mc =
    vtkSmartPointer< vtkMarchingCubes >::New();
#if VTK_MAJOR_VERSION <=5
  mc->SetInput(itk2vtko->GetOutput());
#else
  mc->SetInputData(itk2vtko->GetOutput());
#endif

  mc->SetValue(0, ContourValue);
  mc->Update();

  if (mc->GetOutput()->GetNumberOfCells() == 0)
    {
    std::cerr << "Segmentation failed !" << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "Computing volume of mesh using discretized divergence" << std::endl;
  vtkSmartPointer< vtkMassProperties > mp =
    vtkSmartPointer< vtkMassProperties >::New();
#if VTK_MAJOR_VERSION <=5
  mp->SetInput(mc->GetOutput());
#else
  mp->SetInputData(mc->GetOutput());
#endif
  const double volume = mp->GetVolume();

  std::cout << "Volume of segmentation mm^3 = " << volume << std::endl;

#ifdef LSTK_USE_AIM
  if (!args.GetValueAsString("OutputAIM").empty())
    {
    WriteAIM( args.GetValueAsString("OutputAIM"), ContourValue,
              seg, data, volume );
    }
#endif

  if( args.GetOptionWasSet("OutputMesh"))
    {
    VTK_CREATE(vtkSTLWriter, writer);
    writer->SetFileName(args.GetValueAsString("OutputMesh").c_str());
#if VTK_MAJOR_VERSION <=5
    writer->SetInput( mc->GetOutput() );
#else
    writer->SetInputData( mc->GetOutput() );
#endif
    writer->Write();
    }

  if (args.GetOptionWasSet("Visualize"))
    {
    return ViewImageAndSegmentationSurface(
             image, mc->GetOutput(), roi, args);
    }


  return EXIT_SUCCESS;
}

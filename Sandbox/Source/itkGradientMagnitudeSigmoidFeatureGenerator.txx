/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkGradientMagnitudeSigmoidFeatureGenerator.txx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGradientMagnitudeSigmoidFeatureGenerator_txx
#define __itkGradientMagnitudeSigmoidFeatureGenerator_txx

#include "itkGradientMagnitudeSigmoidFeatureGenerator.h"


namespace itk
{

/*
 * Constructor
 */
template <unsigned int NDimension>
GradientMagnitudeSigmoidFeatureGenerator<NDimension>
::GradientMagnitudeSigmoidFeatureGenerator()
{
  this->SetNumberOfRequiredInputs( 1 );

  this->m_GradientFilter = GradientFilterType::New();
  this->m_SigmoidFilter = SigmoidFilterType::New();

  typename OutputImageSpatialObjectType::Pointer outputObject = OutputImageSpatialObjectType::New();

  this->ProcessObject::SetNthOutput( 0, outputObject.GetPointer() );
}


/*
 * Destructor
 */
template <unsigned int NDimension>
GradientMagnitudeSigmoidFeatureGenerator<NDimension>
::~GradientMagnitudeSigmoidFeatureGenerator()
{
}

template <unsigned int NDimension>
void
GradientMagnitudeSigmoidFeatureGenerator<NDimension>
::SetInput( const SpatialObjectType * spatialObject )
{
  // Process object is not const-correct so the const casting is required.
  this->SetNthInput(0, const_cast<SpatialObjectType *>( spatialObject ));
}

template <unsigned int NDimension>
const typename GradientMagnitudeSigmoidFeatureGenerator<NDimension>::SpatialObjectType *
GradientMagnitudeSigmoidFeatureGenerator<NDimension>
::GetFeature() const
{
  if (this->GetNumberOfOutputs() < 1)
    {
    return 0;
    }

  return static_cast<const SpatialObjectType*>(this->ProcessObject::GetOutput(0));

}


/*
 * PrintSelf
 */
template <unsigned int NDimension>
void
GradientMagnitudeSigmoidFeatureGenerator<NDimension>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
}


/*
 * Generate Data
 */
template <unsigned int NDimension>
void
GradientMagnitudeSigmoidFeatureGenerator<NDimension>
::GenerateData()
{
  typename InputImageSpatialObjectType::ConstPointer inputObject = 
    dynamic_cast<const InputImageSpatialObjectType * >( this->ProcessObject::GetInput(0) );

  if( !inputObject )
    {
    itkExceptionMacro("Missing input spatial object");
    }

  const InputImageType * inputImage = inputObject->GetImage();

  if( !inputImage )
    {
    itkExceptionMacro("Missing input image");
    }

  std::cout << "INPUT IMAGE" << std::endl;
  inputImage->Print( std::cout );

  this->m_GradientFilter->SetInput( inputImage );
  this->m_SigmoidFilter->SetInput( this->m_GradientFilter->GetOutput() );

  this->m_SigmoidFilter->Update();

  typename OutputImageType::Pointer outputImage = this->m_SigmoidFilter->GetOutput();

  outputImage->DisconnectPipeline();

  typedef itk::ImageSpatialObject< Dimension, OutputPixelType > OutputImageSpatialObjectType;

  OutputImageSpatialObjectType * outputObject =
    dynamic_cast< OutputImageSpatialObjectType * >(this->ProcessObject::GetOutput(0));

  outputObject->SetImage( outputImage );
}

} // end namespace itk

#endif
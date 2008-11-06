/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkSinglePhaseLevelSetSegmentationModule.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkSinglePhaseLevelSetSegmentationModule_h
#define __itkSinglePhaseLevelSetSegmentationModule_h

#include "itkSegmentationModule.h"
#include "itkImageSpatialObject.h"

namespace itk
{

/** \class SinglePhaseLevelSetSegmentationModule
 * \brief Class applies a single-phase level set segmentation method
 *
 * SpatialObjects are used as inputs and outputs of this class.
 *
 * \ingroup SpatialObjectFilters
 */
template <unsigned int NDimension>
class ITK_EXPORT SinglePhaseLevelSetSegmentationModule : public SegmentationModule<NDimension>
{
public:
  /** Standard class typedefs. */
  typedef SinglePhaseLevelSetSegmentationModule         Self;
  typedef SegmentationModule<NDimension>                Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SinglePhaseLevelSetSegmentationModule, SegmentationModule);

  /** Dimension of the space */
  itkStaticConstMacro(Dimension, unsigned int, NDimension);

  /** Type of spatialObject that will be passed as input and output of this
   * segmentation method. */
  typedef typename Superclass::SpatialObjectType         SpatialObjectType;
  typedef typename Superclass::SpatialObjectPointer      SpatialObjectPointer;

  /** Types of the input, feature and output images. */
  typedef float                                         InputPixelType;
  typedef float                                         FeaturePixelType;
  typedef float                                         OutputPixelType;
  typedef itk::Image< InputPixelType, NDimension >      InputImageType;
  typedef itk::Image< FeaturePixelType, NDimension >    FeatureImageType;
  typedef itk::Image< OutputPixelType, NDimension >     OutputImageType;

  /** Types of the Spatial objects used for input, feature and output images. */
  typedef itk::ImageSpatialObject< NDimension, InputPixelType >     InputSpatialObjectType;
  typedef itk::ImageSpatialObject< NDimension, FeaturePixelType >   FeatureSpatialObjectType;
  typedef itk::ImageSpatialObject< NDimension, OutputPixelType >    OutputSpatialObjectType;


protected:
  SinglePhaseLevelSetSegmentationModule();
  virtual ~SinglePhaseLevelSetSegmentationModule();
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Method invoked by the pipeline in order to trigger the computation of
   * the segmentation. */
  void  GenerateData ();

private:
  SinglePhaseLevelSetSegmentationModule(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
# include "itkSinglePhaseLevelSetSegmentationModule.txx"
#endif

#endif
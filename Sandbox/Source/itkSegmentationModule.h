/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkSegmentationModule.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkSegmentationModule_h
#define __itkSegmentationModule_h

#include "itkProcessObject.h"
#include "itkImage.h"
#include "itkDataObjectDecorator.h"
#include "itkSpatialObject.h"

namespace itk
{

/** \class SegmentationModule
 * \brief Class that generates features (typically images) used as input for a segmentation method.
 *
 * The typical use of this class would be to generate the edge-map needed by a
 * Level Set filter to internally compute its speed image.
 *
 * SpatialObjects are used as inputs and outputs of this class.
 *
 * \ingroup SpatialObjectFilters
 */
template <unsigned int NDimension>
class ITK_EXPORT SegmentationModule : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef SegmentationModule              Self;
  typedef ProcessObject                 Superclass;
  typedef SmartPointer<Self>            Pointer;
  typedef SmartPointer<const Self>      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SegmentationModule, ProcessObject);

  /** Dimension of the space */
  itkStaticConstMacro(Dimension, unsigned int, NDimension);

  /** Type of spatialObject that will be passed as input and output of this
   * segmentation method. */
  typedef SpatialObject< NDimension >           SpatialObjectType;
  typedef typename SpatialObjectType::Pointer   SpatialObjectPointer;

  /** Input data that will be used for generating the feature. */
  void SetInput( const SpatialObjectType * input );
  const SpatialObjectType * GetInput() const;

  /** Output data that carries the feature in the form of a
   * SpatialObject. */
  const SpatialObjectType * GetFeature() const;


protected:
  SegmentationModule();
  virtual ~SegmentationModule();
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Method invoked by the pipeline in order to trigger the computation of
   * the segmentation. */
  void  GenerateData ();

private:
  SegmentationModule(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
# include "itkSegmentationModule.txx"
#endif

#endif
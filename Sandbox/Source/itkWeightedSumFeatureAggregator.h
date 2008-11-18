/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkWeightedSumFeatureAggregator.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkWeightedSumFeatureAggregator_h
#define __itkWeightedSumFeatureAggregator_h

#include "itkFeatureAggregator.h"

namespace itk
{

/** \class WeightedSumFeatureAggregator
 * \brief Class for combining multiple features into a single one by computing
 * the pixel-wise minimum. 
 *
 * This class generates a new feature object containing an image that is
 * computed as the pixel-wise minimum of all the input feature images.
 *
 * \warning This class assumes that all the images have the same: origin,
 * spacing, orientation, and that they are represented in the same image grid.
 * mixing strategies.
 *
 * SpatialObjects are used as inputs and outputs of this class.
 *
 * \ingroup SpatialObjectFilters
 */
template <unsigned int NDimension>
class ITK_EXPORT WeightedSumFeatureAggregator : public FeatureAggregator<NDimension>
{
public:
  /** Standard class typedefs. */
  typedef WeightedSumFeatureAggregator            Self;
  typedef FeatureAggregator<NDimension>       Superclass;
  typedef SmartPointer<Self>                  Pointer;
  typedef SmartPointer<const Self>            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(WeightedSumFeatureAggregator, FeatureAggregator);

  /** Dimension of the space */
  itkStaticConstMacro(Dimension, unsigned int, NDimension);


protected:
  WeightedSumFeatureAggregator();
  virtual ~WeightedSumFeatureAggregator();
  void PrintSelf(std::ostream& os, Indent indent) const;


private:
  WeightedSumFeatureAggregator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  void ConsolidateFeatures();

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
# include "itkWeightedSumFeatureAggregator.txx"
#endif

#endif
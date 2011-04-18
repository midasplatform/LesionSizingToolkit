/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    itkVectorLevelSetFunction.txx
Language:  C++
Date:      $Date$
Version:   $Revision$

Copyright (c) Insight Software Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVectorLevelSetFunction_txx
#define __itkVectorLevelSetFunction_txx

#include "itkVectorLevelSetFunction.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"

namespace itk {

template <class TImageType>
VectorLevelSetFunction<TImageType>
::VectorLevelSetFunction()
{
  m_EpsilonMagnitude = 1.0e-5;
  m_UseMinimalCurvature = false;
}

template <class TImageType>
typename VectorLevelSetFunction<TImageType>::ScalarValueType
VectorLevelSetFunction<TImageType>
::ComputeCurvatureTerms( const NeighborhoodType &neighborhood,
                         const FloatOffsetType &offset,
                         unsigned int phase,
                         GlobalDataStruct *gd )
{
  ScalarValueType curvature = NumericTraits<ScalarValueType>::Zero;
  for (unsigned int i = 0; i < this->m_NumberOfPhases; i++)
    {
    // Curvature weight due to interaction of phases: 'phase' and 'i'
    ScalarValueType cw = this->m_CurvatureWeights( phase, i );
    if (cw != NumericTraits<ScalarValueType>::Zero)
      {
      curvature += cw * this->ComputeCurvatureTerm( neighborhood, offset,
                                                    i, gd );
      }
    }

  // the cuvature speed can be used to spatially modify the effects of
  // cuvature. FIXME check if this needs more weighting.
  curvature *= this->CurvatureSpeed( neighborhood, offset, phase );

  PhaseDataStruct *pd = &(gd->m_PhaseData[phase]);

  // FIXME valgrind reports an uninitialized var on this line ?
  pd->m_MaxCurvatureChange = vnl_math_max( pd->m_MaxCurvatureChange,
                                           vnl_math_abs(curvature) );

  return curvature;
}

template <class TImageType>
typename VectorLevelSetFunction<TImageType>::ScalarValueType
VectorLevelSetFunction<TImageType>
::ComputeCurvatureTerm( const NeighborhoodType &neighborhood,
                        const FloatOffsetType &offset,
                        unsigned int phase,
                        GlobalDataStruct *gd)
{
  if ( m_UseMinimalCurvature == false )
    {
    return this->ComputeMeanCurvature(neighborhood, offset, phase, gd);
    }
  else
    {
    if (ImageDimension == 3)
      {
      return this->ComputeMinimalCurvature(neighborhood, offset, phase, gd);
      }
    else if (ImageDimension == 2)
      {
      return this->ComputeMeanCurvature(neighborhood, offset, phase, gd);
      }
    else
      {
      return this->ComputeMinimalCurvature(neighborhood, offset, phase, gd);
      }
    }
}


template< class TImageType>
typename VectorLevelSetFunction< TImageType >::ScalarValueType
VectorLevelSetFunction< TImageType >
::ComputeMinimalCurvature(
  const NeighborhoodType &itkNotUsed(neighborhood),
  const FloatOffsetType& itkNotUsed(offset),
  unsigned int phase,
  GlobalDataStruct *gd)
{
  PhaseDataStruct *pd = &(gd->m_PhaseData[phase]);

  unsigned int i, j, n;
  ScalarValueType gradMag = vcl_sqrt(pd->m_GradMagSqr);
  ScalarValueType Pgrad[ImageDimension][ImageDimension];
  ScalarValueType tmp_matrix[ImageDimension][ImageDimension];
  const ScalarValueType ZERO = NumericTraits<ScalarValueType>::Zero;
  vnl_matrix_fixed<ScalarValueType, ImageDimension, ImageDimension> Curve;
  const ScalarValueType MIN_EIG = NumericTraits<ScalarValueType>::min();

  ScalarValueType mincurve;

  for (i = 0; i < ImageDimension; i++)
    {
    Pgrad[i][i] = 1.0 - pd->m_dx[i] * pd->m_dx[i]/gradMag;
    for (j = i+1; j < ImageDimension; j++)
      {
      Pgrad[i][j]= pd->m_dx[i] * pd->m_dx[j]/gradMag;
      Pgrad[j][i] = Pgrad[i][j];
      }
    }

  //Compute Pgrad * Hessian * Pgrad
  for (i = 0; i < ImageDimension; i++)
    {
    for (j = i; j < ImageDimension; j++)
      {
      tmp_matrix[i][j]= ZERO;
      for (n = 0; n < ImageDimension; n++)
        {
        tmp_matrix[i][j] += Pgrad[i][n] * pd->m_dxy[n][j];
        }
      tmp_matrix[j][i]=tmp_matrix[i][j];
      }
    }

  for (i = 0; i < ImageDimension; i++)
    {
    for (j = i; j < ImageDimension; j++)
      {
      Curve(i,j) = ZERO;
      for (n = 0; n < ImageDimension; n++)
        {
        Curve(i,j) += tmp_matrix[i][n] * Pgrad[n][j];
        }
      Curve(j,i) = Curve(i,j);
      }
    }

  //Eigensystem
  vnl_symmetric_eigensystem<ScalarValueType>  eig(Curve);

  mincurve=vnl_math_abs(eig.get_eigenvalue(ImageDimension-1));
  for (i = 0; i < ImageDimension; i++)
    {
    if(vnl_math_abs(eig.get_eigenvalue(i)) < mincurve &&
       vnl_math_abs(eig.get_eigenvalue(i)) > MIN_EIG)
      {
      mincurve = vnl_math_abs(eig.get_eigenvalue(i));
      }
    }

  return ( mincurve / gradMag );
}

template <class TImageType>
typename VectorLevelSetFunction<TImageType>::ScalarValueType
VectorLevelSetFunction<TImageType>::ComputeMeanCurvature(
  const NeighborhoodType &itkNotUsed(neighborhood),
  const FloatOffsetType &itkNotUsed(offset), unsigned int phase, GlobalDataStruct *gd)
{
  PhaseDataStruct *pd = &(gd->m_PhaseData[phase]);

  // Calculate the mean curvature
  ScalarValueType curvature_term = NumericTraits<ScalarValueType>::Zero;
  unsigned int i, j;


  for (i = 0; i < ImageDimension; i++)
    {
    for(j = 0; j < ImageDimension; j++)
      {
      if(j != i)
        {
        curvature_term -= pd->m_dx[i] * pd->m_dx[j] * pd->m_dxy[i][j];
        curvature_term += pd->m_dxy[j][j] * pd->m_dx[i] * pd->m_dx[i];
        }
      }
    }

  return (curvature_term / pd->m_GradMagSqr );
}

template <class TImageType>
typename VectorLevelSetFunction<TImageType>::VectorType
VectorLevelSetFunction<TImageType>::InitializeZeroVectorConstant()
{
  VectorType ans;
  ans.Fill( NumericTraits<ScalarValueType>::Zero );
  return ans;
}

template <class TImageType>
typename VectorLevelSetFunction<TImageType>::VectorType
VectorLevelSetFunction<TImageType>::m_ZeroVectorConstant =
VectorLevelSetFunction<TImageType>::InitializeZeroVectorConstant();

template <class TImageType>
void
VectorLevelSetFunction<TImageType>::
PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "WaveDT: " << m_WaveDT << std::endl;
  os << indent << "DT: " << m_DT << std::endl;
  os << indent << "UseMinimalCurvature " << m_UseMinimalCurvature << std::endl;
  os << indent << "EpsilonMagnitude: " << m_EpsilonMagnitude << std::endl;
  os << indent << "AdvectionWeights: " << m_AdvectionWeights << std::endl;
  os << indent << "PropagationWeights: " << m_PropagationWeights << std::endl;
  os << indent << "CurvatureWeights: " << m_CurvatureWeights << std::endl;
  os << indent << "LaplacianSmoothingWeights: " << m_LaplacianSmoothingWeights << std::endl;
}

template< class TImageType >
double VectorLevelSetFunction<TImageType>::m_WaveDT = 1.0/(2.0 * ImageDimension);

template < class TImageType >
double VectorLevelSetFunction<TImageType>::m_DT     = 1.0/(2.0 * ImageDimension);

template< class TImageType >
typename VectorLevelSetFunction< TImageType >::TimeStepType
VectorLevelSetFunction<TImageType>
::ComputeGlobalTimeStep(void *GlobalData) const
{
  // return the maximum of all the phases.

  TimeStepType dt = 0.0;
  for( unsigned int i = 0; i < this->m_NumberOfPhases; i++ )
    {
    double t = this->ComputeGlobalTimeStep( GlobalData, i );
    dt = (dt > t) ? dt : t;
    }
  return dt;
}

template< class TImageType >
typename VectorLevelSetFunction< TImageType >::TimeStepType
VectorLevelSetFunction<TImageType>
::ComputeGlobalTimeStep(void *GlobalData, unsigned int phase) const
{
  TimeStepType dt;
  GlobalDataStruct *gd = (GlobalDataStruct *)GlobalData;

  PhaseDataStruct *d = &(gd->m_PhaseData[phase]);

  d->m_MaxAdvectionChange += d->m_MaxPropagationChange;

  if (vnl_math_abs(d->m_MaxCurvatureChange) > 0.0)
    {
    if (d->m_MaxAdvectionChange > 0.0)
      {
      dt = vnl_math_min((m_WaveDT / d->m_MaxAdvectionChange),
                        (    m_DT / d->m_MaxCurvatureChange ));
      }
    else
      {
      dt = m_DT / d->m_MaxCurvatureChange;
      }
    }
  else
    {
    if (d->m_MaxAdvectionChange > 0.0)
      {
      dt = m_WaveDT / d->m_MaxAdvectionChange;
      }
    else
      {
      dt = 0.0;
      }
    }

  double maxScaleCoefficient = 0.0;
  for (unsigned int i=0; i<ImageDimension; i++)
    {
    maxScaleCoefficient = vnl_math_max(this->m_ScaleCoefficients[i],maxScaleCoefficient);
    }
  dt /= maxScaleCoefficient;

  // reset the values
  d->m_MaxAdvectionChange   = NumericTraits<ScalarValueType>::Zero;
  d->m_MaxPropagationChange = NumericTraits<ScalarValueType>::Zero;
  d->m_MaxCurvatureChange   = NumericTraits<ScalarValueType>::Zero;

  return dt;
}

template< class TImageType >
void
VectorLevelSetFunction< TImageType>
::Initialize(const RadiusType &r)
{
  this->SetRadius(r);

  // Dummy neighborhood.
  NeighborhoodType it;
  it.SetRadius( r );

  // Find the center index of the neighborhood.
  m_Center =  it.Size() / 2;

  // Get the stride length for each axis.
  for(unsigned int i = 0; i < ImageDimension; i++)
    {  m_xStride[i] = it.GetStride(i); }
}

template< class TImageType >
typename VectorLevelSetFunction< TImageType >::ScalarValueType
VectorLevelSetFunction< TImageType >
::ComputeUpdate( const NeighborhoodType &it,
                 void *globalData,
                 unsigned int phase,
                 const FloatOffsetType& offset )
{
  unsigned int i;
  const PixelType center_value  = it.GetCenterPixel();

  m_NeighborhoodScales = this->ComputeNeighborhoodScales();

  ScalarValueType laplacian_term;
  ScalarValueType propagation_term;
  ScalarValueType curvature_term;
  ScalarValueType advection_term;

  VectorType advection_field;

  // Global data structure
  GlobalDataStruct *gd = (GlobalDataStruct *)globalData;

  for (i = 0; i < this->m_NumberOfPhases; i++)
    {
    gd->m_PhaseData[i].m_GradMagSqr = 1.0e-6;
    }

  // Compute derivative terms etc for the current phase and phases that
  // have crosstalk with it.
  //
  // FIXME Currently we are computing this for every phase. In practice, they
  // need to be computed only for the current phase and for every phase that
  // has crosstalk with the current phase. So a check for zero crosstalk
  // terms needs to be done as a speedup step.
  for (i = 0; i < this->m_NumberOfPhases; i++)
    {
    this->ComputeDerivativesForPhase( it, offset, i, gd );
    }

  curvature_term   = this->ComputeCurvatureTerms  ( it, offset, phase, gd );
  advection_term   = this->ComputeAdvectionTerms  ( it, offset, phase, gd );
  propagation_term = this->ComputePropagationTerms( it, offset, phase, gd );
  laplacian_term   = this->ComputeLaplacianTerms  ( it, offset, phase, gd );

  // Return the combination of all the terms.
  return ( ScalarValueType ) ( curvature_term - propagation_term - advection_term - laplacian_term );
}

template< class TImageType >
void *
VectorLevelSetFunction< TImageType >
::GetGlobalDataPointer() const
{
  if (this->m_NumberOfPhases == 0)
    {
    itkExceptionMacro( <<
      "Number of phases in a VectorLevelSet must be at least 1." );
    }
  GlobalDataStruct *ans = new GlobalDataStruct();

  // Initialize all the phases.
  ans->m_PhaseData = new PhaseDataStruct[this->m_NumberOfPhases];
  for (unsigned int i = 0; i < this->m_NumberOfPhases; i++)
    {
    ans->m_PhaseData[i].m_MaxAdvectionChange   = NumericTraits<ScalarValueType>::Zero;
    ans->m_PhaseData[i].m_MaxPropagationChange = NumericTraits<ScalarValueType>::Zero;
    ans->m_PhaseData[i].m_MaxCurvatureChange   = NumericTraits<ScalarValueType>::Zero;
    }
  return ans;
}


template< class TImageType >
void
VectorLevelSetFunction< TImageType >
::ComputeDerivativesForPhase( const NeighborhoodType &it,
                              const FloatOffsetType &offset,
                              unsigned int pid,
                              GlobalDataStruct *gd )
{
  // FIXME Could this could be computed outside this loop ?
  const PixelType center_value  = it.GetCenterPixel();

  // Compute the Hessian matrix and various other derivatives.  Some of these
  // derivatives may be used by overloaded virtual functions.
  PhaseDataStruct *pd = &(gd->m_PhaseData[pid]);

  for (unsigned int i = 0; i < ImageDimension; i++)
    {
    const unsigned int positionA =
      static_cast<unsigned int>( m_Center + m_xStride[i]);
    const unsigned int positionB =
      static_cast<unsigned int>( m_Center - m_xStride[i]);

    PixelType pixelA = it.GetPixel( positionA );
    PixelType pixelB = it.GetPixel( positionB );

    ScalarValueType pixelAc = pixelA[pid];
    ScalarValueType pixelBc = pixelB[pid];
    ScalarValueType centerVc = center_value[pid];

    pd->m_dx[i] = 0.5 * ( pixelAc - pixelBc ) * m_NeighborhoodScales[i];
    pd->m_dxy[i][i] =   ( pixelAc + pixelBc - 2.0 * centerVc ) * vnl_math_sqr(m_NeighborhoodScales[i]);

    pd->m_dx_forward[i]  = ( pixelAc - centerVc ) * m_NeighborhoodScales[i];
    pd->m_dx_backward[i] = ( centerVc - pixelBc ) * m_NeighborhoodScales[i];
    pd->m_GradMagSqr += pd->m_dx[i] * pd->m_dx[i];

    for( unsigned int j = i+1; j < ImageDimension; j++ )
      {
      const unsigned int positionAa = static_cast<unsigned int>(
        m_Center - m_xStride[i] - m_xStride[j] );
      const unsigned int positionBa = static_cast<unsigned int>(
        m_Center - m_xStride[i] + m_xStride[j] );
      const unsigned int positionCa = static_cast<unsigned int>(
        m_Center + m_xStride[i] - m_xStride[j] );
      const unsigned int positionDa = static_cast<unsigned int>(
        m_Center + m_xStride[i] + m_xStride[j] );

      pd->m_dxy[i][j] = pd->m_dxy[j][i] = 0.25 * ( it.GetPixel( positionAa )[pid]
                                                 - it.GetPixel( positionBa )[pid]
                                                 - it.GetPixel( positionCa )[pid]
                                                 + it.GetPixel( positionDa )[pid] )
                                          * m_NeighborhoodScales[i] * m_NeighborhoodScales[j];
      }
    }
}

template< class TImageType >
void
VectorLevelSetFunction< TImageType >
::ReleaseGlobalDataPointer( void *GlobalData ) const
{
  GlobalDataStruct * gd = (GlobalDataStruct *)GlobalData;
  delete [] gd->m_PhaseData;
  delete gd;
}

template< class TImageType >
typename VectorLevelSetFunction<TImageType>::ScalarValueType
VectorLevelSetFunction< TImageType >
::ComputeAdvectionTerms( const NeighborhoodType &it,
                         const FloatOffsetType &offset,
                         unsigned int phase,
                         GlobalDataStruct * gd) const
{
  const ScalarValueType ZERO = NumericTraits<ScalarValueType>::Zero;
  ScalarValueType advection_term = ZERO;

  for (unsigned int component = 0; 
      component < this->m_NumberOfComponents; component++)
    {
    ScalarValueType w = this->m_AdvectionWeights( phase, component );
    if (w != ZERO)
      {
      VectorType advection_field = this->AdvectionField( it, offset,
                                              component, gd );

      for (unsigned int i = 0; i < ImageDimension; i++)
        {
        ScalarValueType x_energy = w * advection_field[i];

        if (x_energy > ZERO)
          {
          advection_term += advection_field[i] * gd->m_PhaseData[phase].m_dx_backward[i];
          }
        else
          {
          advection_term += advection_field[i] * gd->m_PhaseData[phase].m_dx_forward[i];
          }

        gd->m_PhaseData[phase].m_MaxAdvectionChange
          = vnl_math_max( gd->m_PhaseData[phase].m_MaxAdvectionChange,
                          vnl_math_abs(x_energy));  
        }
      }
    }

  return advection_term;
}

template< class TImageType >
typename VectorLevelSetFunction<TImageType>::ScalarValueType
VectorLevelSetFunction< TImageType >
::ComputePropagationTerms( const NeighborhoodType &it,
                         const FloatOffsetType &offset,
                         unsigned int phase,
                         GlobalDataStruct * gd) const
{
  const ScalarValueType ZERO = NumericTraits<ScalarValueType>::Zero;
  ScalarValueType propagation_term = ZERO;

  for (unsigned int component = 0; 
      component < this->m_NumberOfComponents; component++)
    {

    // Propagation weight due to interaction of phases: 'phase' and 'component'
    ScalarValueType w = this->m_PropagationWeights( phase, component );

    if (w != ZERO)
      {
      ScalarValueType propagationSpeed =
        this->PropagationSpeed(it, offset, component, gd);
      propagation_term += w * propagationSpeed;
      }
    }

  // FIXME  Should we put this in the component loop or outside it ?
  if (propagation_term != ZERO)
    {

    // Construct upwind gradient values for use in the propagation speed term:
    //  $\beta G(\mathbf{x})\mid\nabla\phi\mid$
    //
    // The following scheme for ``upwinding'' in the normal direction is taken
    // from Sethian, Ch. 6 as referenced above.
    //
    ScalarValueType propagation_gradient = ZERO;

    if ( propagation_term > ZERO )
      {
      for (unsigned int i = 0; i< ImageDimension; i++)
        {
        propagation_gradient +=
          vnl_math_sqr( vnl_math_max( gd->m_PhaseData[phase].m_dx_backward[i], ZERO) )
          + vnl_math_sqr( vnl_math_min(gd->m_PhaseData[phase].m_dx_forward[i], ZERO) );
        }
      }
    else
      {
      for (unsigned int i = 0; i< ImageDimension; i++)
        {
        propagation_gradient +=
          vnl_math_sqr( vnl_math_min(gd->m_PhaseData[phase].m_dx_backward[i], ZERO) )
          + vnl_math_sqr( vnl_math_max(gd->m_PhaseData[phase].m_dx_forward[i],  ZERO) );
        }
      }

    // Collect energy change from propagation term.  This will be used in
    // calculating the maximum time step that can be taken for this iteration.
    gd->m_PhaseData[phase].m_MaxPropagationChange =
      vnl_math_max(gd->m_PhaseData[phase].m_MaxPropagationChange,
                   vnl_math_abs(propagation_term));

    propagation_term *= vcl_sqrt( propagation_gradient );
    }

  return propagation_term;
}

template< class TImageType >
typename VectorLevelSetFunction<TImageType>::ScalarValueType
VectorLevelSetFunction< TImageType >
::ComputeLaplacianTerms( const NeighborhoodType &it,
                         const FloatOffsetType &offset,
                         unsigned int phase,
                         GlobalDataStruct * gd) const
{
  const ScalarValueType ZERO = NumericTraits<ScalarValueType>::Zero;
  ScalarValueType laplacian_term = ZERO;

  for (unsigned int component = 0; 
      component < this->m_NumberOfComponents; component++)
    {

    // Propagation weight due to interaction of phases: 'phase' and 'component'
    ScalarValueType w = this->m_LaplacianSmoothingWeights( phase, component );

    if (w != ZERO)
      {
      ScalarValueType laplacian_term_for_component = 
        this->LaplacianSmoothingSpeed( it, offset, component, gd );

      laplacian_term += w * laplacian_term_for_component;
      }
    }

  ScalarValueType laplacian = ZERO;

  // Compute the laplacian for this phase using the existing second 
  // derivative values
  for (unsigned int i = 0;i < ImageDimension; i++)
    {
    laplacian += gd->m_PhaseData[phase].m_dxy[i][i];
    }

  // Scale the laplacian by its speed and weight
  laplacian_term *= laplacian;

  return laplacian_term;
}

} // end namespace itk

#endif
/*
 * Copyright (c) 2023:
 * Author: Dmitry Logashenko, Shuai Lu
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

/*
 * FV discretization for a general conservation law.
 */

/* UG4 headers: */
#include "common/util/provider.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/conv_shape.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{
namespace Conservation_Law_FV{
	
////////////////////////////////////////////////////////////////////////////////
//	check the grid and the shape functions
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
void ConservationLawFV<TDomain>::prepare_setting
(
	const std::vector<LFEID> & vLfeID,
	bool bNonRegular
)
{
//	check the grid
	if (bNonRegular)
		UG_THROW ("ERROR in ConservationLawFV:"
				" Currently, hanging nodes are not supported.\n");

//	check number of the components
	if (vLfeID.size () != 1)
		UG_THROW ("Wrong number of functions: The ElemDisc 'ConservationLawFV'"
					" discretizes a scalar conservation law. Specify only 1 function.");

//	check that these are the Nedelec elements
	if (vLfeID[0].order () != 1 || vLfeID[0].type () != LFEID::LAGRANGE)
		UG_THROW ("Error: ConservationLawFV works with the Lagrange-1 functions only.");
}		

////////////////////////////////////////////////////////////////////////////////
// assembling
////////////////////////////////////////////////////////////////////////////////

/// prepares the loop over the elements
template<typename TDomain>
template<typename TElem>
void ConservationLawFV<TDomain>::prepare_element_loop
(
	ReferenceObjectID roid, ///< only elements with this roid are looped over
	int si ///< and only in this subdomain
)
{
	if(m_spConvShape.invalid())
		UG_THROW("ConservationLawFV::prepare_element_loop:"
						" Upwind has not been set.");
						
	typedef FV1Geometry<TElem, dim> TFVGeom;
	TFVGeom& geo = GeomProvider<TFVGeom>::get ();
	static const int refDim = TElem::dim;
	
	const MathVector<refDim>* vSCVFip = geo.scvf_local_ips();
	const size_t numSCVFip = geo.num_scvf_ips();
	const MathVector<refDim>* vSCVip = geo.scv_local_ips();
	const size_t numSCVip = geo.num_scv_ips();
	m_imDiffusion.template 		set_local_ips<refDim>(vSCVFip,numSCVFip, false);
	m_imFlux.template 			set_local_ips<refDim>(vSCVFip,numSCVFip, false);		
	m_imSource.template 		set_local_ips<refDim>(vSCVip,numSCVip, false);
		//	init upwind for element type
	if(!m_spConvShape->template set_geometry_type<TFVGeom>(geo))
		UG_THROW("ConservationLawFV::prep_elem_loop:"
					" Cannot init upwind for element type.");
	
}

/// finalizes the loop over the elements
template<typename TDomain>
template<typename TElem>
void ConservationLawFV<TDomain>::finish_element_loop ()
{
}

/// prepares a given element for assembling
template<typename TDomain>
template<typename TElem>
void ConservationLawFV<TDomain>::prepare_element
(
	const LocalVector & u, ///< local solution at the dofs associated with elem
	GridObject * elem, ///< element to prepare
	ReferenceObjectID roid,  // id of reference element used for assembling
	const position_type vCornerCoords [] ///< coordinates of the corners of the element
)
{
	typedef FV1Geometry<TElem, dim> TFVGeom;
	
	TElem * pElem = static_cast<TElem*> (elem);
	
// 	update the FV geometry for this element
	TFVGeom& geo = GeomProvider<TFVGeom>::get ();
	try
	{
		geo.update (elem, vCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("ConservationLawFV: Cannot update the Finite Volume Geometry for an element.");
	
	//	set global positions
	const MathVector<dim>* vSCVFip = geo.scvf_global_ips();
	const size_t numSCVFip = geo.num_scvf_ips();
	const MathVector<dim>* vSCVip = geo.scv_global_ips();
	const size_t numSCVip = geo.num_scv_ips();
	m_imDiffusion.			set_global_ips(vSCVFip, numSCVFip);
	m_imFlux.				set_global_ips(vSCVFip, numSCVFip);	
	m_imSource.				set_global_ips(vSCVip, numSCVip);
}

template <class TVector>
static TVector CalculateCenter(GridObject* o, const TVector* coords)
{
	TVector v;
	VecSet(v, 0);

	size_t numCoords = 0;
	switch(o->base_object_id()){
		case VERTEX: numCoords = 1; break;
		case EDGE: numCoords = static_cast<Edge*>(o)->num_vertices(); break;
		case FACE: numCoords = static_cast<Face*>(o)->num_vertices(); break;
		case VOLUME: numCoords = static_cast<Volume*>(o)->num_vertices(); break;
		default: UG_THROW("Unknown element type."); break;
	}

	for(size_t i = 0; i < numCoords; ++i)
		VecAdd(v, v, coords[i]);

	if(numCoords > 0)
		VecScale(v, v, 1. / (number)numCoords);

	return v;
}



/// assembles the stiffness part of the local defect
template<typename TDomain>
template<typename TElem>
void ConservationLawFV<TDomain>::ass_dA_elem
(
	LocalVector & d, ///< the local defect
	const LocalVector & u, ///< local solution at the dofs associated with elem
	GridObject * elem, ///< element to prepare
	const position_type vCornerCoords [] ///< coordinates of the corners of the element
)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	typedef FV1Geometry<TElem, dim> TFVGeom;
	
	TElem * pElem = static_cast<TElem*> (elem);
	
// 	Get finite volume geometry
	const TFVGeom& geo = GeomProvider<TFVGeom>::get ();

    const size_t numSh = geo.num_sh();
	const size_t numScvf = geo.num_scvf();
	

//	get conv shapes
//	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo, false);

	if(m_imDiffusion.data_given()|| m_imFlux.data_given())
	{
	// 	loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < numScvf; ++ip)
		{
		// 	get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

		/////////////////////////////////////////////////////
		// Diffusive Term
		/////////////////////////////////////////////////////
			if(m_imDiffusion.data_given())
			{
			//	to compute D \nabla c
				MathVector<dim> Dgrad_c, grad_c;

			// 	compute gradient and shape at ip
				VecSet(grad_c, 0.0);
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					VecScaleAppend(grad_c, u(_C_,sh), scvf.global_grad(sh));

			//	scale by diffusion tensor
				MatVecMult(Dgrad_c, m_imDiffusion[ip], grad_c);

			// 	Compute flux
				const number diff_flux = VecDot(Dgrad_c, scvf.normal());

			// 	Add to local defect
				d(_C_, scvf.from()) -= diff_flux;
				d(_C_, scvf.to()  ) += diff_flux;
			}
			
		/////////////////////////////////////////////////////
		// Flux Term
		/////////////////////////////////////////////////////			
			if(m_imFlux.data_given())
			{
			//	sum up flux
				const number flux = VecDot(m_imFlux[ip], scvf.normal());

			//  add to local defect
				d(_C_, scvf.from()) += flux;
				d(_C_, scvf.to()  ) -= flux;
				// For debug
				//flux_Flux = flux;
			}
		}//end loop
	
	}
}

/// assembles the local Jacobian of the stiffness part
template<typename TDomain>
template<typename TElem>
void ConservationLawFV<TDomain>::ass_JA_elem
(
	LocalMatrix & J, ///< the local matrix
	const LocalVector & u, ///< local solution at the dofs associated with elem
	GridObject * elem, ///< element to prepare
	const position_type vCornerCoords [] ///< coordinates of the corners of the element
)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	typedef FV1Geometry<TElem, dim> TFVGeom;
	
	TElem * pElem = static_cast<TElem*> (elem);
	
// 	Get finite volume geometry
	const TFVGeom& geo = GeomProvider<TFVGeom>::get ();

    const size_t numSh = geo.num_sh();
	const size_t numScvf = geo.num_scvf();
	
//	Diff. Tensor times Gradient
	MathVector<dim> Dgrad;

//	get conv shapes
//	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo, false);

//	Diffusion and Velocity Term
	if(m_imDiffusion.data_given())
	{
	// 	loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < numScvf; ++ip)
		{
		// 	get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);
		////////////////////////////////////////////////////
		// Diffusive Term
		////////////////////////////////////////////////////
			if(m_imDiffusion.data_given())
			{
				#ifdef UG_ENABLE_DEBUG_LOGS
				//	DID_CONV_DIFF_MP
				number D_diff_flux_sum = 0.0;
				#endif

			// 	loop shape functions
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				{
				// 	Compute Diffusion Tensor times Gradient
					MatVecMult(Dgrad, m_imDiffusion[ip], scvf.global_grad(sh));

				//	Compute flux at IP
					const number D_diff_flux = VecDot(Dgrad, scvf.normal());
					//UG_DLOG(DID_CONV_DIFF_MP, 2, ">>OCT_DISC_DEBUG: " << "local_clfv_impl.h: " << "ass_JA_elem(): " << "sh # "  << sh << " ; normalSize scvf # " << ip << ": " << VecLength(scvf.normal()) << "; \t from "<< scvf.from() << "; to " << scvf.to() << "; D_diff_flux: " << D_diff_flux << "; scvf.global_grad(sh): " << scvf.global_grad(sh) << std::endl);

				// 	Add flux term to local matrix // HIER MATRIXINDIZES!!!
					UG_ASSERT((scvf.from() < J.num_row_dof(_C_)) && (scvf.to() < J.num_col_dof(_C_)),
							  "Bad local dof-index on element with object-id " << elem->base_object_id()
							  << " with center: " << CalculateCenter(elem, vCornerCoords));

					J(_C_, scvf.from(), _C_, sh) -= D_diff_flux;
					J(_C_, scvf.to()  , _C_, sh) += D_diff_flux;

					#ifdef UG_ENABLE_DEBUG_LOGS
					//	DID_CONV_DIFF_MP
					D_diff_flux_sum += D_diff_flux;
					#endif
				}

				//UG_DLOG(DID_CONV_DIFF_MP, 2, "D_diff_flux_sum = " << D_diff_flux_sum << std::endl << std::endl);
			}
		}// end loop
	}
}

/// assembles the mass part of the local defect
template<typename TDomain>
template<typename TElem>
void ConservationLawFV<TDomain>::ass_dM_elem
(
	LocalVector & d, ///< the local defect
	const LocalVector & u, ///< local solution at the dofs associated with elem
	GridObject * elem, ///< element to prepare
	const position_type vCornerCoords [] ///< coordinates of the corners of the element
)
{
}

/// assembles the local Jacobian of the mass part
template<typename TDomain>
template<typename TElem>
void ConservationLawFV<TDomain>::ass_JM_elem
(
	LocalMatrix & J, ///< the local matrix
	const LocalVector & u, ///< local solution at the dofs associated with elem
	GridObject * elem, ///< element to prepare
	const position_type vCornerCoords [] ///< coordinates of the corners of the element
)
{
}

/// assembles the local right-hand side
template<typename TDomain>
template<typename TElem>
void ConservationLawFV<TDomain>::ass_rhs_elem
(
	LocalVector & b,
	GridObject * elem,
	const position_type vCornerCoords []
)
{
	typedef FV1Geometry<TElem, dim> TFVGeom;
	
	TElem * pElem = static_cast<TElem*> (elem);
	
// 	Get finite volume geometry
	const TFVGeom& geo = GeomProvider<TFVGeom>::get ();

	// loop Sub Control Volumes (SCV)
	if ( m_imSource.data_given() ) {
		for ( size_t ip = 0; ip < geo.num_scv(); ++ip ) {
			// get current SCV
			const typename TFVGeom::SCV& scv = geo.scv( ip );

			// get associated node
			const int co = scv.node_id();

			// Add to local rhs
			b(_C_, co) += m_imSource[ip] * scv.volume();
			//UG_LOG("b(_C_, co) = " << b(_C_, co) << "; \t ip " << ip << "; \t co " << co << "; \t scv_vol " << scv.volume() << "; \t m_imSource[ip] " << m_imSource[ip] << std::endl);
		}
	}
}



//	computes the linearized defect w.r.t to the diffusion
template<typename TDomain>
template <typename TElem>
void ConservationLawFV<TDomain>::
lin_def_diffusion(const LocalVector& u,
                  std::vector<std::vector<MathMatrix<dim,dim> > > vvvLinDef[],
                  const size_t nip)
{
	typedef FV1Geometry<TElem, dim> TFVGeom;
	
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo, true);

//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

//  loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

	// 	compute gradient at ip
		MathVector<dim> grad_u;	VecSet(grad_u, 0.0);
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			VecScaleAppend(grad_u, u(_C_,sh), scvf.global_grad(sh));

	//	compute the lin defect at this ip
		MathMatrix<dim,dim> linDefect;

	//	part coming from -\nabla u * \vec{n}
		for(size_t k=0; k < (size_t)dim; ++k)
			for(size_t j = 0; j < (size_t)dim; ++j)
				linDefect(j,k) = (scvf.normal())[j] * grad_u[k];

	//	add contribution from convection shapes
		if(convShape.non_zero_deriv_diffusion())
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				MatAdd(linDefect, convShape.D_diffusion(ip, sh), u(_C_, sh));

	//	add contributions
		vvvLinDef[ip][_C_][scvf.from()] -= linDefect;
		vvvLinDef[ip][_C_][scvf.to()  ] += linDefect;
	}
}

template<typename TDomain>
template <typename TElem>
void ConservationLawFV<TDomain>::
lin_def_flux(const LocalVector& u,
             std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
             const size_t nip)
{
	typedef FV1Geometry<TElem, dim> TFVGeom;
	
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

//  loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

	//	add parts for both sides of scvf
		vvvLinDef[ip][_C_][scvf.from()] += scvf.normal();
		vvvLinDef[ip][_C_][scvf.to()] -= scvf.normal();
	}
}

//	computes the linearized defect w.r.t to the source
template<typename TDomain>
template <typename TElem>
void ConservationLawFV<TDomain>::
lin_def_source(const LocalVector& u,
               std::vector<std::vector<number> > vvvLinDef[],
               const size_t nip)
{
	typedef FV1Geometry<TElem, dim> TFVGeom;
	
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	set lin defect
		vvvLinDef[ip][_C_][co] = scv.volume();
	}
}


////////////////////////////////////////////////////////////////////////////////
//	upwind
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
void ConservationLawFV<TDomain>::
set_upwind(SmartPtr<IConvectionShapes<dim> > shapes) {m_spConvShape = shapes;}


////////////////////////////////////////////////////////////////////////////////
//	Value and Gradient
////////////////////////////////////////////////////////////////////////////////
template<typename TDomain>
template <typename TElem>
void ConservationLawFV<TDomain>::
ex_value(number vValue[],
         const MathVector<dim> vGlobIP[],
         number time, int si,
         const LocalVector& u,
         GridObject* elem,
         const MathVector<dim> vCornerCoords[],
         const MathVector<FV1Geometry<TElem, dim>::dim> vLocIP[],
         const size_t nip,
         bool bDeriv,
         std::vector<std::vector<number> > vvvDeriv[])
{
	typedef FV1Geometry<TElem, dim> TFVGeom;
	
	TElem * pElem = static_cast<TElem*> (elem);
	
// 	update the FV geometry for this element
	TFVGeom& geo = GeomProvider<TFVGeom>::get ();

//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

//	number of shape functions
	static const size_t numSH =	ref_elem_type::numCorners;


//	SCVF ip
	if(vLocIP == geo.scvf_local_ips())
	{
	//	Loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	Get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

		//	compute concentration at ip
			vValue[ip] = 0.0;
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				vValue[ip] += u(_C_, sh) * scvf.shape(sh);

		//	compute derivative w.r.t. to unknowns iff needed
			if(bDeriv)
			{
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					vvvDeriv[ip][_C_][sh] = scvf.shape(sh);

				// do not forget that number of DoFs (== vvvDeriv[ip][_C_])
				// might be > scvf.num_sh() in case of hanging nodes!
				size_t ndof = vvvDeriv[ip][_C_].size();
				for (size_t sh = scvf.num_sh(); sh < ndof; ++sh)
					vvvDeriv[ip][_C_][sh] = 0.0;
			}
		}
	}
//	SCV ip
	else if(vLocIP == geo.scv_local_ips())
	{
	//	Loop Sub Control Volumes (SCV)
		for(size_t ip = 0; ip < geo.num_scv(); ++ip)
		{
		// 	Get current SCV
			const typename TFVGeom::SCV& scv = geo.scv(ip);

		//	get corner of SCV
			const size_t co = scv.node_id();

		//	solution at ip
			vValue[ip] = u(_C_, co);

		//	set derivatives if needed
			if(bDeriv)
			{
				size_t ndof = vvvDeriv[ip][_C_].size();
				for(size_t sh = 0; sh < ndof; ++sh)
					vvvDeriv[ip][_C_][sh] = (sh==co) ? 1.0 : 0.0;
			}
		}
	}
// 	general case
	else
	{
	//	get trial space
		LagrangeP1<ref_elem_type>& rTrialSpace = Provider<LagrangeP1<ref_elem_type> >::get();

	//	storage for shape function at ip
		number vShape[numSH];

	//	loop ips
		for(size_t ip = 0; ip < nip; ++ip)
		{
		//	evaluate at shapes at ip
			rTrialSpace.shapes(vShape, vLocIP[ip]);

		//	compute concentration at ip
			vValue[ip] = 0.0;
			for(size_t sh = 0; sh < numSH; ++sh)
				vValue[ip] += u(_C_, sh) * vShape[sh];

		//	compute derivative w.r.t. to unknowns iff needed
		//	\todo: maybe store shapes directly in vvvDeriv
			if(bDeriv)
			{
				for(size_t sh = 0; sh < numSH; ++sh)
					vvvDeriv[ip][_C_][sh] = vShape[sh];

				// beware of hanging nodes!
				size_t ndof = vvvDeriv[ip][_C_].size();
				for (size_t sh = numSH; sh < ndof; ++sh)
					vvvDeriv[ip][_C_][sh] = 0.0;
			}
		}
	}
}

template<typename TDomain>
template <typename TElem>
void ConservationLawFV<TDomain>::
ex_grad(MathVector<dim> vValue[],
        const MathVector<dim> vGlobIP[],
        number time, int si,
        const LocalVector& u,
        GridObject* elem,
        const MathVector<dim> vCornerCoords[],
        const MathVector<FV1Geometry<TElem, dim>::dim> vLocIP[],
        const size_t nip,
        bool bDeriv,
        std::vector<std::vector<MathVector<dim> > > vvvDeriv[])
{
	typedef FV1Geometry<TElem, dim> TFVGeom;
	
	// 	Get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

//	reference dimension
	static const int refDim = ref_elem_type::dim;

//	number of shape functions
	static const size_t numSH =	ref_elem_type::numCorners;


//	SCVF ip
	if(vLocIP == geo.scvf_local_ips())
	{
	//	Loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	Get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

			VecSet(vValue[ip], 0.0);

			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				VecScaleAppend(vValue[ip], u(_C_, sh), scvf.global_grad(sh));

			if(bDeriv)
			{
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					vvvDeriv[ip][_C_][sh] = scvf.global_grad(sh);

				// beware of hanging nodes!
				size_t ndof = vvvDeriv[ip][_C_].size();
				for (size_t sh = scvf.num_sh(); sh < ndof; ++sh)
					vvvDeriv[ip][_C_][sh] = 0.0;
			}
		}
	}
// 	general case
	else
	{
	//	get trial space
		LagrangeP1<ref_elem_type>& rTrialSpace = Provider<LagrangeP1<ref_elem_type> >::get();

	//	storage for shape function at ip
		MathVector<refDim> vLocGrad[numSH];
		MathVector<refDim> locGrad;

	//	Reference Mapping
		MathMatrix<dim, refDim> JTInv;
		ReferenceMapping<ref_elem_type, dim> mapping(vCornerCoords);

	//	loop ips
		for(size_t ip = 0; ip < nip; ++ip)
		{
		//	evaluate at shapes at ip
			rTrialSpace.grads(vLocGrad, vLocIP[ip]);

		//	compute grad at ip
			VecSet(locGrad, 0.0);
			for(size_t sh = 0; sh < numSH; ++sh)
				VecScaleAppend(locGrad, u(_C_, sh), vLocGrad[sh]);

		//	compute global grad
			mapping.jacobian_transposed_inverse(JTInv, vLocIP[ip]);
			MatVecMult(vValue[ip], JTInv, locGrad);

		//	compute derivative w.r.t. to unknowns iff needed
			if(bDeriv)
			{
				for(size_t sh = 0; sh < numSH; ++sh)
					MatVecMult(vvvDeriv[ip][_C_][sh], JTInv, vLocGrad[sh]);

				// beware of hanging nodes!
				size_t ndof = vvvDeriv[ip][_C_].size();
				for (size_t sh = numSH; sh < ndof; ++sh)
					vvvDeriv[ip][_C_][sh] = 0.0;
			}
		}
	}
};


////////////////////////////////////////////////////////////////////////////////
//	register assembling functions
////////////////////////////////////////////////////////////////////////////////

/// registers the local assembler functions for all the elements and dimensions
template<typename TDomain>
void ConservationLawFV<TDomain>::register_all_loc_discr_funcs ()
{
//	get all grid element types in this dimension and below
	typedef typename domain_traits<dim>::DimElemList ElemList;

//	switch assemble functions
	boost::mpl::for_each<ElemList> (RegisterLocalDiscr (this));
}

/// registers the local assembler functions for a given element
template<typename TDomain>
template<typename TElem> // the element to register for
void ConservationLawFV<TDomain>::register_loc_discr_func ()
{
	static const ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;
	static const int refDim = reference_element_traits<TElem>::dim;
	
	this->clear_add_fct(id);
	
	this->set_prep_elem_loop_fct(id, &T::template prepare_element_loop<TElem>);
	this->set_fsh_elem_loop_fct	(id, &T::template finish_element_loop<TElem>);
	this->set_prep_elem_fct		(id, &T::template prepare_element<TElem>);
	this->set_add_def_A_elem_fct(id, &T::template ass_dA_elem<TElem>);
	this->set_add_jac_A_elem_fct(id, &T::template ass_JA_elem<TElem>);
	this->set_add_def_M_elem_fct(id, &T::template ass_dM_elem<TElem>);
	this->set_add_jac_M_elem_fct(id, &T::template ass_JM_elem<TElem>);
	this->set_add_rhs_elem_fct	(id, &T::template ass_rhs_elem<TElem>);
	
	m_imDiffusion.set_fct(id, this, &T::template lin_def_diffusion<TElem>);
	m_imFlux.set_fct(id, this, &T::template lin_def_flux<TElem>);
	m_imSource.	  set_fct(id, this, &T::template lin_def_source<TElem>);
	
	//	exports
	m_exValue->	   template set_fct<T,refDim>(id, this, &T::template ex_value<TElem>);
	m_exGrad->template set_fct<T,refDim>(id, this, &T::template ex_grad<TElem>);
}

////////////////////////////////////////////////////////////////////////////////
///	class constructor
////////////////////////////////////////////////////////////////////////////////
template<typename TDomain>
ConservationLawFV<TDomain>::
ConservationLawFV
(
	const char * functions, ///< name of the unknown u
	const char * subsets ///< subsets where to assemble
)
:	IElemDisc<TDomain> (functions, subsets),
	m_exValue(new DataExport<number, dim>(functions)),
	m_exGrad(new DataExport<MathVector<dim>, dim>(functions)),
	m_spConvShape(new ConvectionShapesNoUpwind<dim>)
{
//	check number of functions
	if (this->num_fct () != 1)
		UG_THROW ("Wrong number of functions: The ElemDisc 'ConservationLawFV'"
					" discretizes a scalar conservation law. Specify only 1 function.");

// init all imports
	init_imports();

//	default value for mass scale
//	set_mass_scale(1.0);
	
//	register assemble functions
	register_all_loc_discr_funcs ();
}

	
////////////////////////////////////////////////////////////////////////////////
//	user data
////////////////////////////////////////////////////////////////////////////////
template <typename TDomain>
typename ConservationLawFV<TDomain>::NumberExport
ConservationLawFV<TDomain>::
value() {return m_exValue;}

template <typename TDomain>
typename ConservationLawFV<TDomain>::GradExport
ConservationLawFV<TDomain>::
gradient() {return m_exGrad;}

template<typename TDomain>
void ConservationLawFV<TDomain>::
init_imports()
{
	//	register imports
		this->register_import(m_imDiffusion);
		this->register_import(m_imFlux);
		this->register_import(m_imSource);
	
		m_imSource.set_rhs_part();
}



//////// Diffusion

template<typename TDomain>
void ConservationLawFV<TDomain>::
set_diffusion(SmartPtr<CplUserData<MathMatrix<dim, dim>, dim> > user)
{
	m_imDiffusion.set_data(user);
}

template<typename TDomain>
void ConservationLawFV<TDomain>::set_diffusion(number val)
{
	if(val == 0.0) set_diffusion(SmartPtr<CplUserData<MathMatrix<dim, dim>, dim> >());
	else set_diffusion(make_sp(new ConstUserMatrix<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConservationLawFV<TDomain>::set_diffusion(const char* fctName)
{
	set_diffusion(LuaUserDataFactory<MathMatrix<dim,dim>, dim>::create(fctName));
}
template<typename TDomain>
void ConservationLawFV<TDomain>::set_diffusion(LuaFunctionHandle fct)
{
	set_diffusion(make_sp(new LuaUserData<MathMatrix<dim,dim>, dim>(fct)));
}
#endif	


//////// Flux

template<typename TDomain>
void ConservationLawFV<TDomain>::
set_flux(SmartPtr<CplUserData<MathVector<dim>, dim> > user)
{
	m_imFlux.set_data(user);
}

template<typename TDomain>
void ConservationLawFV<TDomain>::set_flux(const std::vector<number>& vVel)
{
	bool bZero = true;
	for(size_t i = 0; i < vVel.size(); ++i){
		if(vVel[i] != 0.0) bZero = false;
	}

	if(bZero) set_flux(SmartPtr<CplUserData<MathVector<dim>, dim> >());
	else set_flux(SmartPtr<ConstUserVector<dim> >(new ConstUserVector<dim>(vVel)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConservationLawFV<TDomain>::
set_flux(const char* fctName)
{
	set_flux(LuaUserDataFactory<MathVector<dim>,dim>::create(fctName));
}
template<typename TDomain>
void ConservationLawFV<TDomain>::
set_flux(LuaFunctionHandle fct)
{
	set_flux(make_sp(new LuaUserData<MathVector<dim>,dim>(fct)));
}
#endif


//////// Source

template<typename TDomain>
void ConservationLawFV<TDomain>::
set_source(SmartPtr<CplUserData<number, dim> > user)
{
	m_imSource.set_data(user);
}

template<typename TDomain>
void ConservationLawFV<TDomain>::
set_source(number val)
{
	if(val == 0.0) set_source(SmartPtr<CplUserData<number, dim> >());
	else set_source(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConservationLawFV<TDomain>::
set_source(const char* fctName)
{
	set_source(LuaUserDataFactory<number,dim>::create(fctName));
}

template<typename TDomain>
void ConservationLawFV<TDomain>::
set_source(LuaFunctionHandle fct)
{
	set_source(make_sp(new LuaUserData<number,dim>(fct)));
}
#endif


////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class ConservationLawFV<Domain1d>;
#endif
#ifdef UG_DIM_2
template class ConservationLawFV<Domain2d>;
#endif
#ifdef UG_DIM_3
template class ConservationLawFV<Domain3d>;
#endif
	
	
} // end namespace Conservation_Law_FV
} // end namespace ug

/* End of File */

/*
 * Copyright (c) 2023:
 * Author: Dmitry Logashenko, Shuai Lu
 * 
 * Based on ConvectionDiffusion
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

	m_imFlux.template 			set_local_ips<refDim>(vSCVFip,numSCVFip, false);		
	m_imSource.template 		set_local_ips<refDim>(vSCVip,numSCVip, false);
	
	m_imVelocity.template 		set_local_ips<refDim>(vSCVFip,numSCVFip, false);
	m_imRelativeK.template 		set_local_ips<refDim>(vSCVip,numSCVip, false);
	
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
	
	m_imFlux.				set_global_ips(vSCVFip, numSCVFip);	
	m_imSource.				set_global_ips(vSCVip, numSCVip);
	m_imVelocity.			set_global_ips(vSCVFip, numSCVFip);	
	m_imRelativeK.			set_global_ips(vSCVip, numSCVip);
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
	
// 	Get finite volume geometry
	const TFVGeom& geo = GeomProvider<TFVGeom>::get ();

    const size_t numSh = geo.num_sh();
	const size_t numScvf = geo.num_scvf();

	// 	loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < numScvf; ++ip)
	{
		// 	get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

		//	sum up flux
		const number flux = VecDot(m_imFlux[ip], scvf.normal());

		//  add to local defect
		d(_C_, scvf.from()) += flux;
		d(_C_, scvf.to()  ) -= flux;
	}//end loop
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
	
// 	Get finite volume geometry
	const TFVGeom& geo = GeomProvider<TFVGeom>::get ();

    const size_t numSh = geo.num_sh();
	const size_t numScvf = geo.num_scvf();


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
const typename ConservationLawFV<TDomain>::conv_shape_type&
ConservationLawFV<TDomain>::
get_updated_conv_shapes(const FVGeometryBase& geo, bool compute_deriv)
{
//	compute upwind shapes for transport equation
//	\todo: we should move this computation into the preparation part of the
//			disc, to only compute the shapes once, reusing them several times.
	if(m_imVelocity.data_given())
	{
	//	get diffusion at ips
		const MathMatrix<dim, dim>* vDiffusion = NULL;
		//if(m_imDiffusion.data_given()) vDiffusion = m_imDiffusion.values();

	//	update convection shapes
		if(!m_spConvShape->update(&geo, m_imVelocity.values(), vDiffusion, compute_deriv))
		{
			UG_LOG("ERROR in 'ConvectionDiffusionMP::get_updated_conv_shapes': "
					"Cannot compute convection shapes.\n");
		}
	}

//	return a const (!!) reference to the upwind
	return *const_cast<const IConvectionShapes<dim>*>(m_spConvShape.get());
}



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

};

////////////////////////////////////////////////////////////////////////////////
//	Export Convection
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
template <typename TElem>
void ConservationLawFV<TDomain>::
ex_conv(MathVector<dim> vValue[],
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
	
	//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo, false);

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
			//const typename TFVGeom::SCV& scv = geo.scv(ip);

			VecSet(vValue[ip], 0.0);

			// 	loop shape functions to get krw at the upwinding node
			number krw = 0.0;
			number factor = 0.0;
			for(size_t sh = 0; sh < convShape.num_sh(); ++sh)
			{
				krw += m_imRelativeK[sh] * convShape(ip, sh);
				factor += convShape(ip, sh);
			}
			if (factor != 0)
				krw = krw / factor;
			//  compute upwind convection
			VecScaleAppend(vValue[ip], krw, m_imVelocity[ip]);

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
	
	m_imFlux.set_fct(id, this, &T::template lin_def_flux<TElem>);
	m_imSource.	  set_fct(id, this, &T::template lin_def_source<TElem>);
	
	//	exports
	m_exValue->template set_fct<T,refDim>(id, this, &T::template ex_value<TElem>);
	m_exGrad->template set_fct<T,refDim>(id, this, &T::template ex_grad<TElem>);
	m_exConv->template set_fct<T,refDim>(id, this, &T::template ex_conv<TElem>);
	
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
	m_exConv(new DataExport<MathVector<dim>, dim>(functions)),
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

template<typename TDomain>
void ConservationLawFV<TDomain>::
init_imports()
{
	//	register imports
	this->register_import(m_imFlux);
	this->register_import(m_imSource);
	this->register_import(m_imVelocity);
	this->register_import(m_imRelativeK);
	
	m_imSource.set_rhs_part();
}


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

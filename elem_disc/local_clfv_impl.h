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
	typedef FV1Geometry<TElem, dim> TFVGeom;

	TFVGeom& geo = GeomProvider<TFVGeom>::get ();
	static const int refDim = TElem::dim;
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
}

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
	
	this->clear_add_fct(id);
	
	this->set_prep_elem_loop_fct(id, & this_type::template prepare_element_loop<TElem>);
	this->set_fsh_elem_loop_fct	(id, & this_type::template finish_element_loop<TElem>);
	this->set_prep_elem_fct		(id, & this_type::template prepare_element<TElem>);
	this->set_add_def_A_elem_fct(id, & this_type::template ass_dA_elem<TElem>);
	this->set_add_jac_A_elem_fct(id, & this_type::template ass_JA_elem<TElem>);
	this->set_add_def_M_elem_fct(id, & this_type::template ass_dM_elem<TElem>);
	this->set_add_jac_M_elem_fct(id, & this_type::template ass_JM_elem<TElem>);
	this->set_add_rhs_elem_fct	(id, & this_type::template ass_rhs_elem<TElem>);
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
:	IElemDisc<TDomain> (functions, subsets)
{
//	check number of functions
	if (this->num_fct () != 1)
		UG_THROW ("Wrong number of functions: The ElemDisc 'ConservationLawFV'"
					" discretizes a scalar conservation law. Specify only 1 function.");

//	register assemble functions
	register_all_loc_discr_funcs ();
}
	
} // end namespace Conservation_Law_FV
} // end namespace ug

/* End of File */

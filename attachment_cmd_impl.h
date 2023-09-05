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

#include <vector>

#include "lib_grid/grid_objects/grid_dim_traits.h"
#include "lib_grid/global_attachments.h"


namespace ug{
namespace Conservation_Law_FV{

/// Copies (scalar) values attached to the full-dimensional elements to a cell-centered grid function.
template <typename TGridFunc>
void CopyGlobAttachmentToGF_CellNumber
(
	const char * attachment_name, ///< name of the global attachment
	TGridFunc * pGF, ///< the grid function to copy to
	size_t fct ///< index of the function in the grid function
)
{	
	typedef typename TGridFunc::domain_type domain_type;
	static const int dim = domain_type::dim;
	typedef typename domain_traits<dim>::grid_base_object elem_type;
	typedef ANumber att_type;
	typedef Grid::AttachmentAccessor<elem_type, att_type> accessor_type;
	
	if (pGF->local_finite_element_id (fct).type () != LFEID::PIECEWISE_CONSTANT)
		UG_THROW ("CopyGlobAttachmentToGF_CellNumber: Not a piecewise-constant grid function specified.");
	
	// get the grid
	auto & grid = * (pGF->domain()->grid ());
	
	// get the attachment
	if (! GlobalAttachments::is_declared (attachment_name))
		UG_THROW ("CopyGlobAttachmentToGF_CellNumber: No global attachment '" << attachment_name << "' found.");
	att_type att = GlobalAttachments::attachment<ANumber> (attachment_name);
	if (! grid.template has_attachment<elem_type> (att))
		UG_THROW ("CopyGlobAttachmentToGF_CellNumber: Global attachment '" << attachment_name << "' does not store element values.");
	accessor_type acc (grid, att);
	
	// loop the elements
	std::vector<DoFIndex> ind (1);
	auto end_iter = pGF->template end<elem_type> ();
	for (auto iter = pGF->template begin<elem_type> (); iter != end_iter; ++iter)
	{
		elem_type * pElem = *iter;
		if (pGF->dof_indices (pElem, fct, ind) != 1)
			UG_THROW ("CopyGlobAttachmentToGF_CellNumber: Non-scalar function specified.");
		DoFRef (*pGF, ind[0]) = acc [pElem];
	}
}

/**
 * Copies (scalar) values attached to the full-dimensional elements to a cell-centered grid function.
 */
template <typename TGridFunc>
void CopyGlobAttachmentToGF_CellNumber
(
	const char * attachment_name, ///< name of the global attachment
	SmartPtr<TGridFunc> spGF, ///< the grid function to copy to
	const char * fct_name ///< name of the function in the grid function
)
{
	FunctionGroup fctGrp;
	try
	{
		fctGrp = spGF->fct_grp_by_name (fct_name);
	}
	UG_CATCH_THROW ("CopyGlobAttachmentToGF_CellNumber: Functions '" << fct_name <<
		"' not all contained in the element approximation space.");
	if (fctGrp.size () != 1)
		UG_THROW ("CopyGlobAttachmentToGF_CellNumber: Only one function component per call supported");
	
	CopyGlobAttachmentToGF_CellNumber (attachment_name, spGF.get (), fctGrp[0]);
}


/// Copies (scalar) values attached to the full-dimensional elements to a nodal grid function.
template <typename TGridFunc>
void CopyGlobAttachmentToGF_NodeNumber
(
	const char * attachment_name, ///< name of the global attachment
	TGridFunc * pGF, ///< the grid function to copy to
	size_t fct ///< index of the function in the grid function
)
{	
	typedef typename TGridFunc::domain_type domain_type;
	static const int dim = domain_type::dim;
	typedef typename geometry_traits<Vertex>::grid_base_object Vertex_type;
	typedef ANumber att_type;
	typedef Grid::AttachmentAccessor<Vertex_type, att_type> accessor_type;
	
	if (pGF->local_finite_element_id (fct).type () != LFEID::LAGRANGE)
		UG_THROW ("CopyGlobAttachmentToGF_NodeNumber: Not a lagrange grid function specified.");
	
	// get the grid
	auto & grid = * (pGF->domain()->grid ());
	
	// get the attachment
	if (! GlobalAttachments::is_declared (attachment_name))
		UG_THROW ("CopyGlobAttachmentToGF_NodeNumber: No global attachment '" << attachment_name << "' found.");
	att_type att = GlobalAttachments::attachment<ANumber> (attachment_name);
	if (! grid.template has_attachment<Vertex_type> (att))
		UG_THROW ("CopyGlobAttachmentToGF_NodeNumber: Global attachment '" << attachment_name << "' does not store element values.");
	accessor_type acc (grid, att);
	
	// loop the nodes
	std::vector<DoFIndex> ind (1);
	auto end_iter = pGF->template end<Vertex_type> ();
	for (auto iter = pGF->template begin<Vertex_type> (); iter != end_iter; ++iter)
	{
		Vertex_type * pNode = *iter;
		if (pGF->dof_indices (pNode, fct, ind) != 1)
			UG_THROW ("CopyGlobAttachmentToGF_NodeNumber: Non-scalar function specified.");
		DoFRef (*pGF, ind[0]) = acc [pNode];
	}
}

/**
 * Copies (scalar) values attached to the full-dimensional elements to a cell-centered grid function.
 */
template <typename TGridFunc>
void CopyGlobAttachmentToGF_NodeNumber
(
	const char * attachment_name, ///< name of the global attachment
	SmartPtr<TGridFunc> spGF, ///< the grid function to copy to
	const char * fct_name ///< name of the function in the grid function
)
{
	FunctionGroup fctGrp;
	try
	{
		fctGrp = spGF->fct_grp_by_name (fct_name);
	}
	UG_CATCH_THROW ("CopyGlobAttachmentToGF_NodeNumber: Functions '" << fct_name <<
		"' not all contained in the element approximation space.");
	if (fctGrp.size () != 1)
		UG_THROW ("CopyGlobAttachmentToGF_NodeNumber: Only one function component per call supported");
	
	CopyGlobAttachmentToGF_NodeNumber (attachment_name, spGF.get (), fctGrp[0]);
}

} // end namespace Conservation_Law_FV
} // end namespace ug

/* End of File */

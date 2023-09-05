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

#ifndef __H__UG__PLUGINS__CONSERVATION_LAW_FV_ATTACHMENT_CMD__
#define __H__UG__PLUGINS__CONSERVATION_LAW_FV_ATTACHMENT_CMD__

#include "common/common.h"
#include "lib_disc/function_spaces/grid_function.h"

namespace ug{
namespace Conservation_Law_FV{

/**
 * Copies (scalar) values attached to the full-dimensional elements to a cell-centered grid function.
 */
template <typename TGridFunc>
void CopyGlobAttachmentToGF_CellNumber
(
	const char * attachment_name, ///< name of the global attachment
	TGridFunc * pGF, ///< the grid function to copy to
	size_t fct ///< index of the function in the grid function
);

/**
 * Copies (scalar) values attached to the full-dimensional elements to a cell-centered grid function.
 */
template <typename TGridFunc>
void CopyGlobAttachmentToGF_CellNumber
(
	const char * attachment_name, ///< name of the global attachment
	SmartPtr<TGridFunc> pGF, ///< the grid function to copy to
	const char * fct_name ///< name of the function in the grid function
);


template <typename TGridFunc>
void CopyGlobAttachmentToGF_NodeNumber
(
	const char * attachment_name, ///< name of the global attachment
	TGridFunc * pGF, ///< the grid function to copy to
	size_t fct ///< index of the function in the grid function
);

/**
 * Copies (scalar) values attached to the full-dimensional elements to a nodal grid function.
 */
template <typename TGridFunc>
void CopyGlobAttachmentToGF_NodeNumber
(
	const char * attachment_name, ///< name of the global attachment
	SmartPtr<TGridFunc> pGF, ///< the grid function to copy to
	const char * fct_name ///< name of the function in the grid function
);

} // end namespace Conservation_Law_FV
} // end namespace ug

#include "attachment_cmd_impl.h"

#endif /* __H__UG__PLUGINS__CONSERVATION_LAW_FV_ATTACHMENT_CMD__ */

/* End of File */

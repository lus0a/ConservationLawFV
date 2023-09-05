/*
 * Copyright (c) 2023
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

/**
 * Plugin of the vertex-centered FV-discretizations of general conservation laws.
 */

/* ug headers: */
#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"

/* discretizations' headers: */
#include "elem_disc/local_clfv.h"

/* auxiliary stuff */
#include "attachment_cmd.h"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace Conservation_Law_FV{

struct Functionality
{
	/**
	 * Function called for the registration of Domain dependent parts.
	 * All Functions and Classes depending on the Domain
	 * are to be placed here when registering. The method is called for all
	 * available Domain types, based on the current build options.
	 *
	 * @param reg	registry
	 * @param grp	group for sorting of functionality
	 */
	template <typename TDomain>
	static void Domain (Registry& reg, string grp)
	{
		string suffix = GetDomainSuffix<TDomain> ();
		string tag = GetDomainTag<TDomain> ();
		static const int dim = TDomain::dim;
		
		{
		typedef ConservationLawFV<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("ConservationLawFV").append (suffix);
		reg.template add_class_<T, TBase > (name, grp)
			.template add_constructor<void (*) (const char*,const char*)> ("Function(s)#Subset(s)")
			.add_method("set_flux", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_flux), "", "Flux")
			#ifdef UG_FOR_LUA
			.add_method("set_flux", static_cast<void (T::*)(const char*)>(&T::set_flux), "", "Flux")
			.add_method("set_flux", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_flux), "", "Flux")
			#endif
			.add_method("set_source", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_source), "", "Source")
			.add_method("set_source", static_cast<void (T::*)(number)>(&T::set_source), "", "Source")
			#ifdef UG_FOR_LUA
			.add_method("set_source", static_cast<void (T::*)(const char*)>(&T::set_source), "", "Source")
			.add_method("set_source", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_source), "", "Source")
			#endif

			.add_method("value", &T::value)
			.add_method("gradient", &T::gradient)

			.set_construct_as_smart_pointer (true);
		reg.add_class_to_group (name, "ConservationLawFV", tag);
		}
	}
	
	/**
	 * Function called for the registration of Domain and Algebra dependent parts.
	 * All Functions and Classes depending on both Domain and Algebra
	 * are to be placed here when registering. The method is called for all
	 * available Domain and Algebra types, based on the current build options.
	 *
	 * @param reg	registry
	 * @param grp	group for sorting of functionality
	 */
	template <typename TDomain, typename TAlgebra>
	static void DomainAlgebra (Registry& reg, string grp)
	{
	//	static const int dim = TDomain::dim;
	//	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra> ();
	//	string tag = GetDomainAlgebraTag<TDomain,TAlgebra> ();
		
		{
		typedef GridFunction<TDomain, TAlgebra> gf_type; 
		reg.add_function("CopyGlobAttachmentToGF_CellNumber", static_cast<void (*)(const char *, SmartPtr<gf_type>, const char *)> (&CopyGlobAttachmentToGF_CellNumber), grp, "Copies attachment to a function", "Attachment#GridFunc#Component");
		}
		
		{
		typedef GridFunction<TDomain, TAlgebra> gf_type; 
		reg.add_function("CopyGlobAttachmentToGF_NodeNumber", static_cast<void (*)(const char *, SmartPtr<gf_type>, const char *)> (&CopyGlobAttachmentToGF_NodeNumber), grp, "Copies attachment to a function", "Attachment#GridFunc#Component");
		}
	};
	
};

} // end namespace Conseration_Law_FV

/**
 * This function is called when the plugin is loaded.
 */
extern "C" void
InitUGPlugin_ConservationLawFV(Registry* reg, string grp)
{
	grp.append("/SpatialDisc/ConservationLawFV");
	typedef Conservation_Law_FV::Functionality Functionality;

	try
	{
		RegisterDomainDependent<Functionality>(*reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // namespace ug

/* End of File */

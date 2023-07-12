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

#ifndef __H__UG__PLUGINS__CONSERVATION_LAW_FV__
#define __H__UG__PLUGINS__CONSERVATION_LAW_FV__

// basic ug4 headers
#include "common/common.h"

// library-specific headers
#include "lib_grid/lg_base.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"
#include "lib_disc/function_spaces/grid_function.h"

/* Discretization's headers: */
#include "lib_disc/spatial_disc/disc_util/conv_shape_interface.h"

namespace ug{
namespace Conservation_Law_FV{

/// FV discretization of the general scalar conservation law
/**
 * This class implements the local FV discretization of the general conservation law:
 * \f{eqnarray*}{
 *  (M (u))_t + \nabla \cdot \mathbf{F} (u) = 0,
 * \f}
 * where
 * <ul>
 * <li> \f$ u \f$		the unknown
 * <li> \f$ M (u) \f$	the mass term (given as a scalar function of u)
 * <li> \f$ \mathbf{F} (u) \f$	the flux function (a vector, possible differential operator of u)
 * </ul>
 *
 * \tparam	TDomain		Domain type
 */
template <typename TDomain>
class ConservationLawFV
	: public IElemDisc<TDomain>
{
private:
///	base class type
	typedef IElemDisc<TDomain> base_type;

///	own type
	typedef ConservationLawFV<TDomain> this_type;

///	domain type
	typedef typename base_type::domain_type domain_type;

///	position type
	typedef typename base_type::position_type position_type;
	
///	world dimension
	static const int dim = base_type::dim;
	
///	abbreviation for the local solution
	static const size_t _C_ = 0;

public:
///	class constructor
	ConservationLawFV
	(
		const char * function, ///< name of the unknown u
		const char * subsets ///< subsets where to assemble
	);
	
///	to compute the upwind shapes
	void set_upwind(SmartPtr<IConvectionShapes<dim> > shapes);
	
protected:
	void init_imports();
	
	/// method to compute the upwind shapes
		SmartPtr<IConvectionShapes<dim> > m_spConvShape;

	///	returns the updated convection shapes
		typedef IConvectionShapes<dim> conv_shape_type;
		const IConvectionShapes<dim>& get_updated_conv_shapes(const FVGeometryBase& geo, bool compute_deriv);
		
			///	computes the concentration
		template <typename TElem>
		void ex_value(number vValue[],
		              const MathVector<dim> vGlobIP[],
		              number time, int si,
		              const LocalVector& u,
		              GridObject* elem,
		              const MathVector<dim> vCornerCoords[],
		              const MathVector<FV1Geometry<TElem, dim>::dim> vLocIP[],
		              const size_t nip,
		              bool bDeriv,
		              std::vector<std::vector<number> > vvvDeriv[]);

	///	computes the gradient of the concentration
		template <typename TElem>
		void ex_grad(MathVector<dim> vValue[],
		             const MathVector<dim> vGlobIP[],
		             number time, int si,
		             const LocalVector& u,
		             GridObject* elem,
		             const MathVector<dim> vCornerCoords[],
		             const MathVector<FV1Geometry<TElem, dim>::dim> vLocIP[],
		             const size_t nip,
		             bool bDeriv,
		             std::vector<std::vector<MathVector<dim> > > vvvDeriv[]);
	
public:
///	sets the diffusion tensor
	/**
	 * This method sets the Diffusion tensor used in computations. If no
	 * Tensor is set, a zero value is assumed.
	 */
	///	\{
		void set_diffusion(SmartPtr<CplUserData<MathMatrix<dim, dim>, dim> > user);
		void set_diffusion(number val);
#ifdef UG_FOR_LUA
		void set_diffusion(const char* fctName);
		void set_diffusion(LuaFunctionHandle fct);
#endif
	///	\}
	
	
	///	sets the flux
	/**
	 * This method sets the Flux. If no field is provided a zero
	 * value is assumed.
	 */
	/// \{
		void set_flux(SmartPtr<CplUserData<MathVector<dim>, dim> > user);
		void set_flux(const std::vector<number>& vVel);
#ifdef UG_FOR_LUA
		void set_flux(const char* fctName);
		void set_flux(LuaFunctionHandle fct);
#endif
	/// \}	
	
	
	///	sets the source / sink term
	/**
	 * This method sets the source/sink value. A zero value is assumed as
	 * default.
	 */
	///	\{
		void set_source(SmartPtr<CplUserData<number, dim> > user);
		void set_source(number val);
#ifdef UG_FOR_LUA
		void set_source(const char* fctName);
		void set_source(LuaFunctionHandle fct);
#endif
	///	\}
	
protected:
	///	Data import for Diffusion
		DataImport<MathMatrix<dim,dim>, dim> m_imDiffusion;
		
	///	Data import for the Flux
		DataImport<MathVector<dim>, dim > m_imFlux;
		
	///	Data import for the right-hand side (volume)
		DataImport<number, dim> m_imSource;

public:
		typedef SmartPtr<CplUserData<number, dim> > NumberExport;
		typedef SmartPtr<CplUserData<MathVector<dim>, dim> > GradExport;

	///	returns the export of the value of associated unknown function
		SmartPtr<CplUserData<number, dim> > value();

	///	returns the export of the gradient of associated unknown function
		SmartPtr<CplUserData<MathVector<dim>, dim> > gradient();

protected:
	///	Export for the concentration
		SmartPtr<DataExport<number, dim> > m_exValue;

	///	Export for the gradient of concentration
		SmartPtr<DataExport<MathVector<dim>, dim> > m_exGrad;
		

//---- Local discretization interface: ----
private:
	
///	check type of the grid and the trial space
	virtual void prepare_setting
	(
		const std::vector<LFEID> & vLfeID,
		bool bNonRegular
	);

//---- Assembling functions: ----
	
	template <typename TElem>
	void prepare_element_loop(ReferenceObjectID roid, int si);

	template <typename TElem>
	void prepare_element(const LocalVector& u, GridObject* elem, ReferenceObjectID roid, const position_type vCornerCoords[]);

	template <typename TElem>
	void finish_element_loop();

	template <typename TElem>
	void ass_JA_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const position_type vCornerCoords[]);

	template <typename TElem>
	void ass_JM_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const position_type vCornerCoords[]);

	template <typename TElem>
	void ass_dA_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const position_type vCornerCoords[]);

	template <typename TElem>
	void ass_dM_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const position_type vCornerCoords[]);

	template <typename TElem>
	void ass_rhs_elem(LocalVector& d, GridObject* elem, const position_type vCornerCoords[]);
	
	
protected:
	///	computes the linearized defect w.r.t to the diffusion
		template <typename TElem>
		void lin_def_diffusion(const LocalVector& u,
		                       std::vector<std::vector<MathMatrix<dim,dim> > > vvvLinDef[],
		                       const size_t nip);

	///	computes the linearized defect w.r.t to the flux
		template <typename TElem>
		void lin_def_flux(const LocalVector& u,
		                  std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
		                  const size_t nip);
						  
	///	computes the linearized defect w.r.t to the source term
		template <typename TElem>
		void lin_def_source(const LocalVector& u,
		                    std::vector<std::vector<number> > vvvLinDef[],
		                    const size_t nip);

//---- Registration of the template functions: ----
private:
	
	void register_all_loc_discr_funcs();

	struct RegisterLocalDiscr {
			RegisterLocalDiscr(this_type * pThis) : m_pThis(pThis){}
			this_type * m_pThis;
			template< typename TElem > void operator() (TElem &)
			{m_pThis->register_loc_discr_func<TElem> ();}
	};

	template <typename TElem>
	void register_loc_discr_func ();


}; // end class ConservationLawFV


} // end namespace Conservation_Law_FV
} // end namespace ug

#include "local_clfv_impl.h"

#endif /* __H__UG__PLUGINS__CONSERVATION_LAW_FV__ */

/* End of File */

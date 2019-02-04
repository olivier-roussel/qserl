#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <qserl/rod3d/rod.h>

#include <eigenpy/eigenpy.hpp>

using namespace boost::python;

typedef Eigen::Matrix<double,6,1> Vector6d;
typedef Eigen::Matrix<double,7,1> Vector7d;

namespace qserl {
  namespace rod3d {
    WorkspaceIntegratedState::IntegrationResultT bind_integrateStateFromBaseWrench
      (RodShPtr rod,
       const Vector6d& v_wrench,
       const Vector7d& v_basePos,
       const WorkspaceIntegratedState::IntegrationOptions& integrationOptions)
      {
        Eigen::Wrenchd wrench (v_wrench);
        //Eigen::Wrenchd wrench;
        //Eigen::Displacementd basePos (v_basePos);
        Eigen::Displacementd basePos; basePos.get() = v_basePos;
        rod->integrateStateFromBaseWrench (wrench, basePos, integrationOptions);
      }

    Vector7d WorkspaceState_node (const WorkspaceState& ws, std::size_t i)
    {
      return ws.nodes()[i].get();
    }

    void exposeToPython()
    {
      class_< std::vector<double> >("StdVector_double")
        .def(vector_indexing_suite<std::vector<double> >());

      enum_ <Parameters::RodModelT> ("RodModelT")
        .value ("RM_INEXTENSIBLE"        , Parameters::RM_INEXTENSIBLE)
        .value ("RM_EXTENSIBLE_SHEARABLE", Parameters::RM_EXTENSIBLE_SHEARABLE)
        .value ("RM_NUMBER_OF_ROD_MODELS", Parameters::RM_NUMBER_OF_ROD_MODELS)
        ;
      enum_ <WorkspaceIntegratedState::IntegrationResultT> ("IntegrationResultT")
        .value ("IR_VALID"                        , WorkspaceIntegratedState::IR_VALID)
        .value ("IR_SINGULAR"                     , WorkspaceIntegratedState::IR_SINGULAR)
        .value ("IR_UNSTABLE"                     , WorkspaceIntegratedState::IR_UNSTABLE)
        .value ("IR_OUT_OF_WRENCH_BOUNDS"         , WorkspaceIntegratedState::IR_OUT_OF_WRENCH_BOUNDS)
        .value ("IR_NUMBER_OF_INTEGRATION_RESULTS", WorkspaceIntegratedState::IR_NUMBER_OF_INTEGRATION_RESULTS)
        ;
      class_<Parameters> ("Parameters", init<>())
        .def_readwrite  ("radius"               ,   &Parameters::radius               )
        .def_readwrite  ("length"               ,   &Parameters::length               )
        // need eigenpy
        .def_readwrite  ("stiffnessCoefficients",   &Parameters::stiffnessCoefficients)
        .def_readwrite  ("rodModel"             ,   &Parameters::rodModel             )
        .def_readwrite  ("numNodes"             ,   &Parameters::numNodes             )
        .def ("setIsotropicStiffnessCoefficientsFromElasticityParameters", &Parameters::setIsotropicStiffnessCoefficientsFromElasticityParameters)
        ;
      class_<Rod, RodShPtr, boost::noncopyable> ("Rod", no_init)
        .def ("create", &Rod::create)
        .staticmethod("create")
        .def ("integrateStateFromBaseWrench", bind_integrateStateFromBaseWrench)
        .def ("integratedState", &Rod::integratedState)
        ;

      class_<WorkspaceState, WorkspaceStateShPtr, boost::noncopyable> ("WorkspaceState", no_init)
        .def ("numNodes", &WorkspaceState::numNodes)
        .def ("node", WorkspaceState_node)
        //.def ("nodes", &WorkspaceState::nodes)
        //.def ("nodesAbsolute6DPositions", &WorkspaceState::nodesAbsolute6DPositions)
        ;
      class_<WorkspaceIntegratedState, WorkspaceIntegratedStateShPtr, bases<WorkspaceState>, boost::noncopyable> ("WorkspaceIntegratedState", no_init)
        .def ("getMMatrix", &WorkspaceIntegratedState::getMMatrix, return_value_policy<return_by_value>())
        .def ("getJMatrix", &WorkspaceIntegratedState::getJMatrix, return_value_policy<return_by_value>())
        .def ("J_det", &WorkspaceIntegratedState::J_det, return_value_policy<return_by_value>())
        .def ("J_nu_sv", &WorkspaceIntegratedState::J_nu_sv, return_value_policy<return_by_value>())
        ;

      class_ <WorkspaceIntegratedState::IntegrationOptions> ("IntegrationOptions", init<>())
        .def_readwrite ("computeJ_nu_sv"  , &WorkspaceIntegratedState::IntegrationOptions::computeJ_nu_sv)
        .def_readwrite ("stop_if_unstable", &WorkspaceIntegratedState::IntegrationOptions::stop_if_unstable)
        .def_readwrite ("keepMuValues"    , &WorkspaceIntegratedState::IntegrationOptions::keepMuValues)
        .def_readwrite ("keepJdet"        , &WorkspaceIntegratedState::IntegrationOptions::keepJdet)
        .def_readwrite ("keepMMatrices"   , &WorkspaceIntegratedState::IntegrationOptions::keepMMatrices)
        .def_readwrite ("keepJMatrices"   , &WorkspaceIntegratedState::IntegrationOptions::keepJMatrices)
        ;
    }
  }
}

namespace eigenpy {
  void exposeMissingMatrices ()
  {
    ENABLE_SPECIFIC_MATRIX_TYPE(Vector6d);
    ENABLE_SPECIFIC_MATRIX_TYPE(Vector7d);

    typedef Eigen::Matrix<double,6,6> Matrix6d;
    ENABLE_SPECIFIC_MATRIX_TYPE(Matrix6d);
  }
}

BOOST_PYTHON_MODULE(rod3d)
{
  boost::python::import ("eigenpy");

  eigenpy::exposeMissingMatrices();

  qserl::rod3d::exposeToPython();
}

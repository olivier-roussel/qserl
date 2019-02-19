#include <boost/python.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <Eigen/Geometry>

#include <qserl/rod3d/rod.h>
#include <qserl/rod3d/ik.h>
#include <qserl/util/explog.h>

#include <eigenpy/eigenpy.hpp>

using namespace boost::python;

using Eigen::Vector3d;
using Eigen::Matrix3d;
typedef Eigen::Matrix<double,6,1> Vector6d;
typedef Eigen::Matrix<double,7,1> Vector7d;
typedef Eigen::Matrix<double,6,6> Matrix6d;

namespace qserl {
  namespace rod3d {
    list displacementToTQ (const Displacement& d)
    {
      list ret;
      // Translation
      for (int i = 0; i < 3; ++i) ret.append (d(i,3));
      // Rotation
      Eigen::Quaterniond q (d.topLeftCorner<3,3>());
      for (int i = 0; i < 4; ++i) ret.append (q.coeffs()[i]);
      return ret;
    }

    Displacement WorkspaceState_node (const WorkspaceState& ws, std::size_t i)
    {
      return ws.nodes()[i];
    }

    Eigen::Vector3d _log3_a(const Eigen::Matrix3d & R)
    {
      return log3 (R);
    }

    tuple _log3_b(const Eigen::Matrix3d & R)
    {
      double t;
      Eigen::Vector3d w = log3 (R, t);
      return make_tuple (w, t);
    }

    Displacement _exp6 (const Vector6d& v) { return exp6 (v); }
    Matrix3d _exp3 (const Vector3d& v) { return exp3 (v); }
    Vector6d _log6 (const Displacement& v) { return log6 (v); }

    void exposeToPython()
    {
      typedef return_value_policy<return_by_value> policy_by_value;

      class_< std::vector<double> >("StdVector_double")
        .def(vector_indexing_suite<std::vector<double> >());
      // Does not work and I don't know why.
      class_< Displacements >("Displacements")
        .def(vector_indexing_suite<Displacements >());

      def ("displacementToTQ", displacementToTQ);
      def ("inv"  , inv<double>);
      def ("exp3" , _exp3);
      def ("log3" , _log3_a);
      def ("log3t", _log3_b);
      def ("exp6" , _exp6);
      def ("log6" , _log6);

      {
        scope parameters =
          class_<Parameters> ("Parameters", init<>())
          .def_readwrite ("radius"               , &Parameters::radius               )
          .def_readwrite ("stiffnessCoefficients", &Parameters::stiffnessCoefficients)
          .def_readwrite ("rodModel"             , &Parameters::rodModel             )
          .def_readwrite ("numNodes"             , &Parameters::numNodes             )
          .def_readwrite ("gravity"              , &Parameters::gravity              )
          .def_readwrite ("unitaryMass"          , &Parameters::unitaryMass          )
          .def_readwrite ("integrationTime"      , &Parameters::integrationTime      )
          .def ("setIsotropicStiffnessCoefficientsFromElasticityParameters", &Parameters::setIsotropicStiffnessCoefficientsFromElasticityParameters)
          ;
        enum_ <Parameters::RodModelT> ("RodModelT") ;
        // Make RodModelT values accessible with Parameters.value
        parameters.attr ("RM_INEXTENSIBLE"              ) = Parameters::RM_INEXTENSIBLE;
        parameters.attr ("RM_EXTENSIBLE_SHEARABLE"      ) = Parameters::RM_EXTENSIBLE_SHEARABLE;
        parameters.attr ("RM_INEXTENSIBLE_WITH_GRAVITY" ) = Parameters::RM_INEXTENSIBLE_WITH_GRAVITY;
        parameters.attr ("RM_NUMBER_OF_ROD_MODELS"      ) = Parameters::RM_NUMBER_OF_ROD_MODELS;
      }

      class_<Rod, RodShPtr, boost::noncopyable> ("Rod", no_init)
        .def ("create", &Rod::create)
        .staticmethod("create")
        .def ("parameters", &Rod::parameters, policy_by_value())
        .def ("integrateStateFromBaseWrench", &Rod::integrateStateFromBaseWrench)
        .def ("integratedState", &Rod::integratedState)
        ;

      class_<WorkspaceState, WorkspaceStateShPtr, boost::noncopyable> ("WorkspaceState", no_init)
        .def ("numNodes", &WorkspaceState::numNodes)
        .def ("node", WorkspaceState_node)
        .def ("nodes", &WorkspaceState::nodes, policy_by_value())
        .def ("base", (const Displacement& (WorkspaceState::*)() const)&WorkspaceState::base, policy_by_value())
        // .def ("nodesAbsolute6DPositions", &WorkspaceState::nodesAbsolute6DPositions)
        ;
      {
        scope wis =
          class_<WorkspaceIntegratedState, WorkspaceIntegratedStateShPtr, bases<WorkspaceState>, boost::noncopyable> ("WorkspaceIntegratedState", no_init)
          .def ("wrench"    , &WorkspaceIntegratedState::wrench)
          .def ("isStable"  , &WorkspaceIntegratedState::isStable)
          .def ("getMMatrix", &WorkspaceIntegratedState::getMMatrix, policy_by_value())
          .def ("getJMatrix", &WorkspaceIntegratedState::getJMatrix, policy_by_value())
          .def ("J_det"     , &WorkspaceIntegratedState::J_det     , policy_by_value())
          .def ("J_nu_sv"   , &WorkspaceIntegratedState::J_nu_sv   , policy_by_value())
          ;
        enum_ <WorkspaceIntegratedState::IntegrationResultT> ("IntegrationResultT");
        // Make IntegrationResultT values accessible with WorkspaceIntegratedState.value
        wis.attr ("IR_VALID"                        ) = WorkspaceIntegratedState::IR_VALID;
        wis.attr ("IR_SINGULAR"                     ) = WorkspaceIntegratedState::IR_SINGULAR;
        wis.attr ("IR_UNSTABLE"                     ) = WorkspaceIntegratedState::IR_UNSTABLE;
        wis.attr ("IR_OUT_OF_WRENCH_BOUNDS"         ) = WorkspaceIntegratedState::IR_OUT_OF_WRENCH_BOUNDS;
        wis.attr ("IR_NUMBER_OF_INTEGRATION_RESULTS") = WorkspaceIntegratedState::IR_NUMBER_OF_INTEGRATION_RESULTS;

        class_ <WorkspaceIntegratedState::IntegrationOptions> ("IntegrationOptions", init<>())
          .def_readwrite ("computeJ_nu_sv"  , &WorkspaceIntegratedState::IntegrationOptions::computeJ_nu_sv)
          .def_readwrite ("stop_if_unstable", &WorkspaceIntegratedState::IntegrationOptions::stop_if_unstable)
          .def_readwrite ("keepMuValues"    , &WorkspaceIntegratedState::IntegrationOptions::keepMuValues)
          .def_readwrite ("keepJdet"        , &WorkspaceIntegratedState::IntegrationOptions::keepJdet)
          .def_readwrite ("keepMMatrices"   , &WorkspaceIntegratedState::IntegrationOptions::keepMMatrices)
          .def_readwrite ("keepJMatrices"   , &WorkspaceIntegratedState::IntegrationOptions::keepJMatrices)
          ;
      }

      {
        scope iks =
          class_ <InverseKinematics> ("InverseKinematics", init<RodShPtr>())
          .def ("compute", &InverseKinematics::compute)
          .add_property ("errorThreshold"  , &InverseKinematics::getErrorThreshold, &InverseKinematics::setErrorThreshold)
          .add_property ("verbosity"       , &InverseKinematics::getVerbosity, &InverseKinematics::setVerbosity)
          .add_property ("maxIterations"   , &InverseKinematics::getMaxIter, &InverseKinematics::setMaxIter)
          .add_property ("scale"           , &InverseKinematics::getScale, &InverseKinematics::setScale)
          .def ("lastIntegrationResult"    , &InverseKinematics::lastIntegrationResult)
          ;
        enum_ <InverseKinematics::ResultT> ("ResultT");
        // Make ResultT values accessible with InverseKinematics.value
        iks.attr ("IK_VALID"                        ) = InverseKinematics::IK_VALID;
        iks.attr ("IK_JACOBIAN_SINGULAR"            ) = InverseKinematics::IK_JACOBIAN_SINGULAR;
        iks.attr ("IK_INTEGRATION_FAILED"           ) = InverseKinematics::IK_INTEGRATION_FAILED;
        iks.attr ("IK_MAX_ITER_REACHED"             ) = InverseKinematics::IK_MAX_ITER_REACHED;
        iks.attr ("IR_NUMBER_OF_INTEGRATION_RESULTS") = InverseKinematics::IR_NUMBER_OF_INTEGRATION_RESULTS;
      }
    }
  }
}

namespace eigenpy {
  void exposeMissingMatrices ()
  {
    typedef qserl::rod3d::Wrench Wrench;
    typedef qserl::rod3d::Displacement Displacement;
    ENABLE_SPECIFIC_MATRIX_TYPE(Wrench);
    ENABLE_SPECIFIC_MATRIX_TYPE(Displacement);

    ENABLE_SPECIFIC_MATRIX_TYPE(Vector7d);

    ENABLE_SPECIFIC_MATRIX_TYPE(Matrix6d);
  }
}

BOOST_PYTHON_MODULE(rod3d)
{
  boost::python::import ("eigenpy");

  eigenpy::exposeMissingMatrices();

  qserl::rod3d::exposeToPython();
}

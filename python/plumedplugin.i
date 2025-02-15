%module openmmplumed

%import(module="openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"
%include "std_string.i"
%include "std_vector.i"

namespace std {
  %template(vectord) vector<double>;
}

%{
#include "PlumedForce.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
#include <numpy/arrayobject.h>
#include <mpi.h>

int isNumpyAvailable() {
    return true;
}
%}

%pythoncode %{
import openmm as mm
%}

%include /home/clay/software/pytorch_build/conda/pytorch/2.3.1/lib/python3.11/site-packages/mpi4py/include/mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);

namespace PlumedPlugin {

class PlumedForce : public OpenMM::Force {
public:
    PlumedForce(const std::string& script);
    PlumedForce(const std::string& script, const MPI_Comm& comm, const int& rank_mod);
    const std::string& getScript() const;
    bool usesPeriodicBoundaryConditions() const;
    void setTemperature(double temperature);
    double getTemperature() const;
    void setMasses(const std::vector<double>& masses);
    const std::vector<double>& getMasses() const;
    void setRestart(bool restart);
    bool getRestart() const;
};

}

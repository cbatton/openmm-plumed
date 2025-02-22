/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2016 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "openmm/OpenMMException.h"
#include "PlumedForce.h"
#include "internal/PlumedForceImpl.h"

using namespace PlumedPlugin;
using namespace OpenMM;
using namespace std;

PlumedForce::PlumedForce(const string& script) : script(script), temperature(-1),
    logStream(stdout), restart(false) {
}

PlumedForce::PlumedForce(const string& script, const MPI_Comm& comm, const int& rank_mod = 1) :
    script(script), temperature(-1), logStream(stdout), restart(false) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    bool in_new_comm = (rank % rank_mod == 0);
    MPI_Comm_split(comm, in_new_comm, rank, &sub_comm);
    MPI_Comm_split(comm, rank, 0, &fake_comm);
    if (in_new_comm) {
        int new_rank, new_size;
        MPI_Comm_rank(sub_comm, &new_rank);
        MPI_Comm_size(sub_comm, &new_size);
    }

}

const string& PlumedForce::getScript() const {
    return script;
}

void PlumedForce::setTemperature(double temperature_) {
    temperature = temperature_;
}

double PlumedForce::getTemperature() const {
    return temperature;
}

void PlumedForce::setMasses(const std::vector<double>& masses_) {
    masses = masses_;
}

const std::vector<double>& PlumedForce::getMasses() const {
    return masses;
}

void PlumedForce::setLogStream(FILE* stream) {

    if (!stream)
        throw OpenMMException("PlumedForce::setLogStream: the stream has to be open");

    logStream = stream;
}

FILE* PlumedForce::getLogStream() const {
    return logStream;
}

ForceImpl* PlumedForce::createImpl() const {
    return new PlumedForceImpl(*this);
}

void PlumedForce::setRestart(bool restart_) {
    restart = restart_;
}

bool PlumedForce::getRestart() const {
    return restart;
}

const MPI_Comm& PlumedForce::getMPIComm() const {
    return sub_comm;
}

const MPI_Comm& PlumedForce::getMPICommFake() const {
    return fake_comm;
}

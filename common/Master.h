/***************************************************************************
 *   Copyright (C) 2009 by Florian Goth   *
 *   fgoth@wthp095   *
 *                                                                         *
 *   All rights reserved.                                                  *
 *                                                                         *
 *   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met: *
 *     * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. *
 *     * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. *
 *     * Neither the name of the <ORGANIZATION> nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission. *
 *                                                                         *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS   *
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT     *
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR *
 *   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR *
 *   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, *
 *   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,   *
 *   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    *
 *   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF *
 *   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING  *
 *   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 *   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.          *
 ***************************************************************************/
#ifndef MASTER_H
#define MASTER_H
#include <sys/time.h>
#include <string>
#include <vector>
#include <list>
#include <inttypes.h>
#include <fstream>
#include "ParallelServer.h"
#include "Models.h"

struct Parameters
{
  typedef float BetaType;
  typedef float UType;
  typedef float VType;
  typedef float EDotType;
    uint32_t N;
    double t;
  UType U;
  BetaType beta;
  double alpha;
    double mu;
    double t_exp;
    double delta;
    
  VType V;
    double W;
  EDotType ed;
    
    bool has_spin_symmetry;
    bool data_is_complex;
    char* data;
    uint datalen;
    
    uint Nx;
    uint Nb;
    double lambda;
    
    uint32_t N_wait;
    uint32_t prebins;
    std::string path;
    std::string binpath;
    std::string idpath;
    double contourlen;
    Models model;
    double delta_s;//the resolution of the time-axis. If t_exp == 0 n_max= beta/delta_s
                   //if t_exp != 0 n_max = t_exp/delta_s
    unsigned int functionpoints;
    bool signiscomplex;
    bool is_Impurity_model;
};

template <class SignType>
class ClientState; //forward declaration

struct Early_ObservableProperties
{
  int idx;
  bool covariance;
  Early_ObservableProperties() : covariance(false) {}
  Early_ObservableProperties(int i) : idx(i), covariance(false) {}
  Early_ObservableProperties(int i, bool c) : idx(i), covariance(c) {}
};

class Master : public ParallelServer
{
public:
    Master(int argc, char *argv[]);
    ~Master();
    template <typename SignType>
    void run(std::vector<Parameters>::const_iterator, std::list<ClientState<SignType> >&);
    std::vector<Parameters> params;
private:
    unsigned int cycles;
    bool usetimer;
    bool analyzeonly;
    struct timeval tstart;
    double interval;
    double elapsedtime();
    void connectclient(const Client& c, const Parameters& curparams);
    void receiveMeasurements();
    template <typename T>
    char* readinDatafromFiles(class RegistryDB&, uint, uint&);
    template <typename SignType>
    void wait_for_all_clients_to_stop(std::list<ClientState<SignType> >& clients, std::ofstream&);
    std::vector<Early_ObservableProperties> observablenames;
    template <typename SignType>
    inline int handleActivity(std::list<ClientState<SignType> >& clients, std::ofstream&);
    template <typename ClientContainer>
    bool wait_for_Init_Ack(const Client& c, ClientContainer& clients, std::ofstream& additionalinfo);
};
#endif

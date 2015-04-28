/***************************************************************************
 *   Copyright (C) 2009 - 2014 by Florian Goth   *
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

#include "ParallelServer.h"
#include <iostream>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <complex>
#include <ctime>
#include <fstream>
#include <valarray>
#include <sys/stat.h>
#include "registry.h"
#include "commands.h"
#include "AverageSign.h"
#include "ClientState.h"

#ifdef _OPENMP
#include <omp.h>
#endif

template <typename T>
struct GetTrait<std::vector<T> >//helper trait to break up a string at various predefined seperators
{
    static std::vector<T> Convert(std::string& arg)
    {
        const string delim("; ,");
        std::vector<T> retval;
        std::size_t pos = 0;
        std::size_t posold = 0;
        do
        {
            retval.push_back(GetTrait<T>::Convert(arg.substr(posold, ( pos = arg.find_first_of(delim, posold) ) - posold) ));
        } while ((posold = arg.find_first_not_of(delim, pos) ) != string::npos);
        return retval;
    }
};

template <>
struct GetTrait<std::vector<std::string> >//helper trait to break up a string at various predefined seperators
{
    static std::vector<std::string> Convert(std::string arg)
    {
        const string delim("; ,");
        std::vector<std::string> retval;
        std::size_t pos = 0;
        std::size_t posold = 0;
        do
        {
            retval.push_back(arg.substr(posold, ( pos = arg.find_first_of(delim, posold) ) - posold) );
        } while ((posold = arg.find_first_not_of(delim, pos) ) != string::npos);
        return retval;
    }
};

enum ObservableDataType {//with examples...
    SCALAR,//<n>
    FUNCTION,//<E>(t)
    VECTORFUNCTION//<n_k>(t)
};

struct ObservableProperties
{
    std::string name;
    bool iscomplex;
    ObservableDataType type;
    bool isk_space;
    bool makes_only_sense_for_impurity;
};

//Note that the names are fixed between Master and Zombie
static const ObservableProperties propertycache[] =
{
    {"AverageOrder", false, SCALAR, false, false},
    {"ParticleNumber", false, SCALAR, false, false},
    {"TotalDoubleOccupancy", false, SCALAR, false, false},
    {"Hybridization_up", true, SCALAR, false, true},
    {"Hybridization_down", true, SCALAR, false, true},
    {"KineticEnergy", false, FUNCTION, false, false},
    {"Magnetization", false, FUNCTION, false, false},
    {"EtaPairing", true, VECTORFUNCTION, true, false},
    {"SpinSpinCorrelation", false, VECTORFUNCTION, false, false},
    {"ChargeChargeCorrelation", false, VECTORFUNCTION, false, false},
    {"kSpaceDensity", true, VECTORFUNCTION, true, false},
    {"LocalDensityVariance", false, VECTORFUNCTION, false, false},
    {"Greensfunction", false, VECTORFUNCTION, false, false},
    {"DensityDensityStructureFactor", true, VECTORFUNCTION, true, false},
    {"DoubleOccupancy", false, VECTORFUNCTION, false, false},
    {"ChargeChargeCorrelatedPart", false, VECTORFUNCTION, false, false},
    {"SpinSpinCorrelatedPart", false, VECTORFUNCTION, false, false},
    {"EtaPairing_Real", false, VECTORFUNCTION, false, false},
    {"Conductance", false, SCALAR, false, true},
    {"SmoothImaginaryGreensfunction", false, VECTORFUNCTION, false, false},
    {"SmoothImaginaryGreensfunction_averaged", false, VECTORFUNCTION, false, false},
    {"MatsubaraFrequencyGreensfunction", true, VECTORFUNCTION, false, false},
    {"MatsubaraFrequencyGreensfunction_averaged", true, VECTORFUNCTION, false, false},
    {"MatsubaraFrequencyGreensfunctionDerivative", true, VECTORFUNCTION, false, false},
    {"SpinSpinCorrelatedPart_Y", false, VECTORFUNCTION, false, false},
    {"SpinSpinCorrelatedPart_X", false, VECTORFUNCTION, false, false},
    {"ImaginarySpinSpinCorrelation_Z", false, VECTORFUNCTION, false, false},
    {"SpinSusceptibility_X", false, SCALAR, false, false},
    {"SpinSusceptibility_Z", false, SCALAR, false, false},
    {"LocalBathGreensfunctions_up", true, VECTORFUNCTION, false, true},
    {"LocalBathGreensfunctions_down", true, VECTORFUNCTION, false, true},
    {"LocalBathGreensfunctions_averaged", true, VECTORFUNCTION, false, true},
    {"KondoCloud_Z", true, FUNCTION, false, true},
    {"KondoCloud_X", true, FUNCTION, false, true},
    {"DiagonalGreensfunction_kspace", true, VECTORFUNCTION, true, false},
    {"OffDiagonalGreensfunction_kspace", true, VECTORFUNCTION, true, false},
    {"SplusSminus_kspace", false, VECTORFUNCTION, true, false},
    {"SpinSpinCorrelation_kspace", false, VECTORFUNCTION, true, false},
    {"ChargeChargeCorrelation_kspace", false, VECTORFUNCTION, true, false},
    {"Gplusplus_kspace", true, VECTORFUNCTION, true, false},
    {"Gminusminus_kspace", false, VECTORFUNCTION, true, false},
};

static const unsigned int nr_of_known_observables = sizeof(propertycache)/sizeof(ObservableProperties);

//Here we store the properties of the model. It's ordered in the same way as the model enums.
static const bool modelhascomplexsign[] =
{
  //Hubbard_chain, cold_atoms, siam, imag_model_from_file, kondo_imp_ti, rashba_chain, exponential Rashba, power law rashba
    false, false, false, true, true, false, true, true
};

using namespace std;
// General tool to strip spaces from both ends:
static inline std::string trim(const std::string& s)
{
    //from "Thinking in C++" by Bruce Eckel
    if (s.length() == 0)
        return s;
    int b = s.find_first_not_of(" \t");
    int e = s.find_last_not_of(" \t");
    if (b == -1) // No non-spaces
        return "";
    return std::string(s, b, e - b + 1);
}

/**
Function object for erasing a client
*/
template <class SignType>
class Client_Eraser
{
public:
    Client_Eraser(const Client& a) : c(a) {}
    inline bool operator()(const ClientState<SignType>& a) const
    {
        return a.client == c;
    }
private:
    const Client& c;
};

void Master::receiveMeasurements()
{}

void Master::connectclient(const Client& c, const Parameters& curparams)
{
    cout<<"A new soul awaits its calling from the master... take this undead seed...."<<endl;
    int32_t seed = lrand48();
    send<char>(c, INIT);//set the zombie's state to init
    send<int32_t>(c, curparams.model);
    send<char>(c, (curparams.t_exp == 0.0 ? 0: 1));//Note to self: we really send it twice...
    switch (curparams.model)
    {
    case HUBBARD_CHAIN:
    case COLD_ATOMS:
        send<uint32_t>(c, curparams.N);//send particle number
        break;
    case SIAM:
        send<double>(c, curparams.V);//send hybridization
        send<double>(c, curparams.W);//send bandwidth
        send<double>(c, curparams.ed);//send dot level
        break;
    case IMAG_MODEL_FROM_FILE:
    {
        send<char>(c, curparams.has_spin_symmetry?1:0);
        send<char>(c, curparams.data_is_complex?1:0);
        send<uint32_t>(c, curparams.N);//send particle number
//	send<int32_t>(c, curparams.datalen);
        std::valarray<char> data(curparams.data, curparams.datalen);
        send<std::valarray<char> >(c, data);
    }
    break;
    case KONDO_IMP_TI:
        send<double>(c, curparams.V);//send hybridization
        send<double>(c, curparams.ed);//send dot level
        send<uint32_t>(c, curparams.Nx);
        send<uint32_t>(c, curparams.Nb);
        send<double>(c, curparams.lambda);//send coupling
        break;
    case RASHBA_CHAIN:
      send<uint32_t>(c, curparams.N);//send particle number
      send<double>(c, curparams.lambda);//send strength of Rashba coupling
      break;
    case RASHBA_CHAIN_EXPONENTIAL:
    case RASHBA_CHAIN_POWER_LAW:
      send<uint32_t>(c, curparams.N);//send particle number
      send<double>(c, curparams.lambda);//send strength of Rashba coupling
      send<double>(c, curparams.alpha);//send strength of decay
      break;
    }
    send<double>(c, curparams.t);//send kinetic energy
    send<double>(c, curparams.U);//send U
    send<double>(c, curparams.beta);
    send<double>(c, curparams.mu);
    send<double>(c, curparams.t_exp);
    send<double>(c, curparams.delta);
    send<uint32_t>(c, curparams.N_wait);//send the waiting time between the measurements
    send<uint32_t>(c, curparams.prebins);//send number of bins to pre-average
    send<double>(c, curparams.delta_s);//send the time resolution
    send<int32_t>(c, seed);//send the seeds
    //now send all the used observables
    send<uint32_t>(c, observablenames.size());//send the number of observables we want to measure
    for (std::vector<Early_ObservableProperties>::const_iterator it = observablenames.begin(); it != observablenames.end(); ++it)//send the observable names
        send<string>(c, propertycache[it->idx].name);
    send<char>(c, START);//finished sending configuration. let the zombie run!
    return;
}

template <typename ClientContainer>
bool Master::wait_for_Init_Ack(const Client& c, ClientContainer& clients, std::ofstream& additionalinfo)
{
    time_t rawtime;
    time ( &rawtime );
    additionalinfo<<"Taking control of: "<<c.getName()<<" at "<<ctime (&rawtime)<<std::endl;
    while(true)
    {
    try
    {
        checkActivity();//see if sth. is happening on any port
    }
    catch (ConnectionError& e)
    {
        const string msg( e.client.getName() + string(" decapitated."));
	time_t rawtime;
        time ( &rawtime );
        cout<<"["<<ctime(&rawtime)<<"] "<<"Activity on Client Socket but no Data received! Assuming Zombie finally dead!"<<endl;
        cout<<msg<<endl;
        additionalinfo<<msg<<std::endl;
        clients.erase(remove_if(clients.begin(), clients.end(), Client_Eraser<typename ClientContainer::value_type::Signtype>(e.client)), clients.end());//remove the client from our array
        removeClient(e.client);//remove the client from the net subsystem
	if(e.client == c) return false;
    }
      if(anyMessagefrom(c))
      {
	if(peeknextSize(c) == 1)
	{
	  Letter<char> ack(recvfromspecificClient(c));
	  if (ack.msg == INIT_RECEIVED) break;
	}
	else
	{
	  discardmessage(c);//can't be an ack. hence we discard that message
	}
      }
    }
    std::cout<<"ACK from "<<c.getName()<<" received!"<<std::endl;
    return true;
}


template <typename T>
void readDatafromFile(T* data, std::string& filename, uint N, uint nr_of_points)
{
    ifstream file(filename.c_str());
    for (uint k = 0; k < N; ++k)
    {
        uint temp;
        file >>std::skipws >> temp;//ignore site number
        for (uint i = 0; i < nr_of_points; ++i)
        {
            double tau;
            file >>std::skipws >>tau;//ignore tau
            file >>std::skipws >> data[k*nr_of_points + i];
        }
    }
}

template <typename T>
char* Master::readinDatafromFiles(RegistryDB& registry, uint N, uint& len)
{
    /**
    Used format in RAM
    spin-sectors, spinsectors consist of sites,
    a single site consists of a bunch of real numbers
    */
    bool spin_symmetry = registry.Get<bool>("modelparameters", "Imag_Model_from_File", "Has_Spin_Symmetry");
    double dtau = registry.Get<double>("modelparameters", "Imag_Model_from_File", "Delta_Tau");
    uint nr_of_tau_points = static_cast<uint>(round(registry.Get<float>("modelparameters" , "Common", "beta")/ dtau)) + 1;
    uint nr_of_points_per_sector = N*nr_of_tau_points;
    len = nr_of_points_per_sector * (spin_symmetry?1:2);
    T* retval = new T[len];
    if (spin_symmetry)
    {
        std::cout<<"Not implemented!"<<std::endl;
        exit(-1);
    }
    else
    {//For now we expect the different spin sectors to reside in two files
        std::string upfile = registry.Get<string>("modelparameters", "Imag_Model_from_File", "Up_Spin_File");
        std::string downfile = registry.Get<string>("modelparameters", "Imag_Model_from_File", "Down_Spin_File");
        readDatafromFile<T>(retval, upfile, N, nr_of_tau_points);
        readDatafromFile<T>(retval + nr_of_points_per_sector, downfile, N, nr_of_tau_points);
    }
    len *= sizeof(T);
    return reinterpret_cast<char*>(retval);
}

double Master::elapsedtime()
{
    timeval curtime;
    gettimeofday(&curtime, NULL);
    return (curtime.tv_sec - tstart.tv_sec) + (curtime.tv_usec - tstart.tv_usec) / 1000000.0;
}

Master::Master(int argc, char *argv[])
try
:
    ParallelServer(argc, argv), cycles(0), usetimer(false), analyzeonly(false), interval(0)
{
    if (argc < 2)
    {
        cout<<"not enough command line arguments!"<<endl;
        exit(-1);
    }
    //decode cmd line
    string toggle(trim(string(argv[argc -2])));
    if (toggle == "-a")
    {
        analyzeonly = true;
    }
    else if (toggle == "-m")
    {
        usetimer = true;//we expect the interval of time to be in minutes
        interval = (atoi(argv[argc-1]))*60.0;//nr of seconds to simulate
        gettimeofday(&tstart, NULL);
    }
    else
    {
        if (toggle == "-n")
            cycles = atoi(argv[argc-1]);
        else
            throw std::runtime_error(string("switch not recognized: ") + toggle );
    }
#ifdef _OPENMP
    cout<<"Nr. of Processors: "<<omp_get_num_procs()<<endl;
#endif
    //initialize the Servers PRNGs with the time
    srand48 (time (0));
    //retrieve data from the config file
    RegistryDB registry("./cfgs");
    //Get model string and decode it to the enum used throughout the programm
    std::string modelstr = registry.Get<std::string>("qmc", "Main", "model");
    std::transform(modelstr.begin(), modelstr.end(), modelstr.begin(), ::toupper);
    Models mymodel;
    if (modelstr == "HUBBARD_CHAIN")
        mymodel = HUBBARD_CHAIN;
    else if (modelstr == "COLD_ATOMS")
        mymodel = COLD_ATOMS;
    else if (modelstr == "SIAM")
        mymodel = SIAM;
    else if (modelstr == "IMAG_MODEL_FROM_FILE")
        mymodel = IMAG_MODEL_FROM_FILE;
    else if (modelstr == "KONDO_IMP_TI")
        mymodel = KONDO_IMP_TI;
    else if (modelstr == "RASHBA_CHAIN")
        mymodel = RASHBA_CHAIN;
    else if (modelstr == "RASHBA_CHAIN_EXPONENTIAL")
        mymodel = RASHBA_CHAIN_EXPONENTIAL;
    else if (modelstr == "RASHBA_CHAIN_POWER_LAW")
        mymodel = RASHBA_CHAIN_POWER_LAW;
    else
    {
        std::cout<<"Model "<<modelstr<<" not recognized. Valid models are: Hubbard_chain, Cold_atoms, Siam, Imag_Model_from_File, Kondo_Imp_TI, Rashba_Chain, Rashba_Chain_exponential, Rashba_Chain_Power_Law"<<std::endl;
        exit(-1);
    }
    //we need to derive the number of different parameter sets that we simulate
    std::vector<Parameters::BetaType> betas = registry.Get<std::vector<float> >("modelparameters", "Common", "beta");
    std::vector<Parameters::UType> us = registry.Get<std::vector<float> >("modelparameters", "Common", "u");
    uint nr_of_sets = betas.size() * us.size();
    std::vector<Parameters::VType> vs;
    std::vector<Parameters::EDotType> eds;
    std::vector<Parameters::BetaType> lambdas;//we declare the vector in case we use it for the Kondo_Imp_TI measurements in the future
    if (mymodel == SIAM)
    {
        vs = registry.Get<std::vector<float> >("modelparameters", "SIAM", "V");
        eds = registry.Get<std::vector<float> >("modelparameters", "SIAM", "e_dot");
        nr_of_sets *= vs.size()*eds.size();
    }
    if (mymodel == KONDO_IMP_TI)
    {
        vs = registry.Get<std::vector<float> >("modelparameters", "Kondo_Imp_TI", "V");
        eds = registry.Get<std::vector<float> >("modelparameters", "Kondo_Imp_TI", "e_dot");
        nr_of_sets *= vs.size()*eds.size();
    }
    if ((mymodel == RASHBA_CHAIN) || (mymodel == RASHBA_CHAIN_EXPONENTIAL) || (mymodel == RASHBA_CHAIN_POWER_LAW))
    {
      lambdas = registry.Get<std::vector<float> >("modelparameters", "RashbaChain", "lambda");
      nr_of_sets *= lambdas.size();
    }
    params.resize(nr_of_sets);
    if ((mymodel == SIAM) || (mymodel == KONDO_IMP_TI))
    {
        for (uint v_idx = 0; v_idx < vs.size(); ++v_idx)
            for (uint ed_idx = 0; ed_idx < eds.size(); ++ed_idx)
                for (uint beta_idx = 0; beta_idx < betas.size(); ++beta_idx)
                    for (uint u_idx = 0; u_idx < us.size(); ++u_idx)
                    {
                        uint cur_idx = v_idx *eds.size() * betas.size()*us.size() + ed_idx*betas.size()*us.size() + beta_idx * us.size() + u_idx;
                        params[cur_idx].beta = betas[beta_idx];
                        params[cur_idx].U = us[u_idx];
                        params[cur_idx].V = vs[v_idx];
                        params[cur_idx].ed = eds[ed_idx];
			params[cur_idx].is_Impurity_model = true;
                    }
    }
    else
    {
      if (mymodel == RASHBA_CHAIN)
      {
	for (uint beta_idx = 0; beta_idx < betas.size(); ++beta_idx)
	  for (uint u_idx = 0; u_idx < us.size(); ++u_idx)
	    for (uint l_idx = 0; l_idx < lambdas.size(); ++l_idx)
	    {
	      uint cur_idx = u_idx * betas.size()*lambdas.size() + beta_idx*lambdas.size() + l_idx;
	      params[cur_idx].beta = betas[beta_idx];
	      params[cur_idx].U = us[u_idx];
	      params[cur_idx].lambda = lambdas[l_idx];
	      params[cur_idx].is_Impurity_model = false;
	    }
      }
      else
      {
        for (uint beta_idx = 0; beta_idx < betas.size(); ++beta_idx)
            for (uint u_idx = 0; u_idx < us.size(); ++u_idx)
            {
	        uint cur_idx = u_idx * betas.size() + beta_idx;
                params[cur_idx].beta = betas[beta_idx];
                params[cur_idx].U = us[u_idx];
		params[cur_idx].is_Impurity_model = false;
            }
      }
    }
    float t_exp = 0.0;
    try
    {
        t_exp = registry.Get<float>("modelparameters", "Common", "t_exp");
    }
    catch (Registry_Block_Data_not_found_Exception&)
        {}//No need to do anything if we don't find the key, we assume that t_exp is zero then,
    for (uint k = 0; k < params.size(); ++k)
    {
        params[k].t_exp = t_exp;
        params[k].model = mymodel;
        params[k].signiscomplex = modelhascomplexsign[mymodel] || (t_exp != 0.);//determine if there is a reason to expect a complex sign
        //retrieve the other parameters
//    curparams.U = registry.Get<float>("modelparameters", "Common", "u");
        params[k].delta = registry.Get<float>("modelparameters", "Common", "delta");
//    curparams.beta = registry.Get<std::vector<float> >("modelparameters", "Common", "beta");
        std::string basepath = registry.Get<string>("modelparameters", "Common", "path");
        params[k].t = registry.Get<float>("modelparameters", "Common", "t");
        params[k].mu = registry.Get<float>("modelparameters", "Common", "my");
        switch (mymodel)
        {
        case HUBBARD_CHAIN:
        case COLD_ATOMS:
            params[k].N = registry.Get<unsigned int>("modelparameters", "HubbardChain", "sites");
            break;
        case SIAM:
            params[k].N = 1;//eases a couple of things...
//        curparams.V = registry.Get<float>("modelparameters", "SIAM", "V");
            params[k].W = registry.Get<float>("modelparameters", "SIAM", "W");
//        curparams.ed = registry.Get<float>("modelparameters", "SIAM", "e_dot");
            break;
        case IMAG_MODEL_FROM_FILE:
            params[k].N = registry.Get<unsigned int>("modelparameters" , "Imag_Model_from_File", "sites");
            params[k].data_is_complex = registry.Get<bool>("modelparameters" , "Imag_Model_from_File", "Data_Is_complex");
            params[k].has_spin_symmetry = registry.Get<bool>("modelparameters" , "Imag_Model_from_File", "Has_Spin_Symmetry");
            if (params[k].data_is_complex)
            {
                params[k].data = readinDatafromFiles<std::complex<double> >(registry, params[k].N, params[k].datalen);
                params[k].signiscomplex = true;
            }
            else
                params[k].data = readinDatafromFiles<double>(registry, params[k].N, params[k].datalen);
            break;
        case KONDO_IMP_TI:
            params[k].N = 1;//eases a couple of things...
//        curparams.V = registry.Get<float>("modelparameters", "Kondo_Imp_TI", "V");
//        curparams.ed = registry.Get<float>("modelparameters", "Kondo_Imp_TI", "e_dot");
            params[k].Nx = registry.Get<unsigned int>("modelparameters", "Kondo_Imp_TI", "Nx");
            params[k].Nb = registry.Get<unsigned int>("modelparameters", "Kondo_Imp_TI", "Nb");
            params[k].lambda = registry.Get<float>("modelparameters", "Kondo_Imp_TI", "lambda");
            break;
	case RASHBA_CHAIN:
	case RASHBA_CHAIN_EXPONENTIAL:
	case RASHBA_CHAIN_POWER_LAW:
	    params[k].N = registry.Get<unsigned int>("modelparameters" , "RashbaChain", "sites");
	    params[k].lambda = registry.Get<double>("modelparameters", "RashbaChain", "lambda");
            break;
        }
        if((mymodel == RASHBA_CHAIN_EXPONENTIAL) || (mymodel == RASHBA_CHAIN_POWER_LAW) )
	{
	  params[k].alpha = registry.Get<double>("modelparameters", "RashbaChain_LR", "alpha");
	}
        params[k].N_wait = registry.Get<unsigned int>("qmc" , "Main", "N_Wait");
	params[k].prebins = registry.Get<unsigned int>("qmc" , "Main", "prebins");
	try
	{
        params[k].delta_s = registry.Get<double>("observables", "Common", "delta_s");
	params[k].functionpoints = static_cast<unsigned int>(floor(params[k].t_exp == 0.?params[k].beta/params[k].delta_s : params[k].t_exp/params[k].delta_s)) + 1;
	}
	catch (Registry_Key_not_found_Exception& e)
        {
        params[k].functionpoints = registry.Get<unsigned int>("observables", "Common", "functionpoints");
	params[k].delta_s = (params[k].t_exp == 0. ? params[k].beta : params[k].t_exp) / params[k].functionpoints;
        }
        params[k].contourlen = 2.0*t_exp + params[k].beta;
        params[k].path = basepath;
        params[k].binpath = registry.Get<string>("qmc" , "Main", "binpath");
        if (*(params[k].binpath.end()-1) != '/')
            params[k].binpath += '/';
//build up the string of the path we use for storing the bins
        params[k].idpath = registry.Get<std::string>("qmc", "Main", "model");
        switch (mymodel)
        {
        case HUBBARD_CHAIN:
        case COLD_ATOMS:
        case IMAG_MODEL_FROM_FILE:
            params[k].idpath += string("N") + toString(params[k].N);
            break;
        case SIAM:
            params[k].idpath += string("_V") + toString(params[k].V) + string("_W") + toString(params[k].W) + string("_ed") + toString(params[k].ed);
            break;
        case KONDO_IMP_TI:
            params[k].idpath += string("_V") + toString(params[k].V) + string("_ed") + toString(params[k].ed);
            params[k].idpath += string("_Nx") + toString(params[k].Nx) + string("_Nb") + toString(params[k].Nb) + string("_lambda") + toString(params[k].lambda);
            break;
	case RASHBA_CHAIN:
	case RASHBA_CHAIN_EXPONENTIAL:
	case RASHBA_CHAIN_POWER_LAW:
	    params[k].idpath += string("_N") + toString(params[k].N) + string("_lambda") + toString(params[k].lambda);
	  break;
        }
        if((mymodel == RASHBA_CHAIN_EXPONENTIAL) || (mymodel == RASHBA_CHAIN_POWER_LAW) )
	  params[k].idpath += string("_alpha") + toString(params[k].alpha);
        params[k].idpath += string("_beta") + toString(params[k].beta) + string("_mu") + toString(params[k].mu) + string("_t_exp") + toString(params[k].t_exp) + string("delta") + toString(params[k].delta) + string("_U") + toString(params[k].U);
        std::string path = basepath + params[k].idpath + string("/");
        //Done reading config files
        mode_t process_mask = umask(0);
        //create output directory
        mkdir(path.c_str(), S_IRWXU | S_IRGRP);
        params[k].binpath = params[k].binpath + params[k].idpath;
        if (*(params[k].binpath.end() - 1) != '/')
            params[k].binpath += '/';
        if (*(params[k].idpath.end() - 1) != '/')
            params[k].idpath += '/';
        //create path for storing the bins
        mkdir(params[k].binpath.c_str(), S_IRWXU | S_IRGRP);
        umask(process_mask);
        int copysuccesful = system((string("cp -r ./cfgs ") + path).c_str());
    }
    bool my_model_is_impurity_model = ((mymodel == SIAM) || (mymodel == KONDO_IMP_TI));
//decode the observables to use
    std::vector<std::string> observable_ids = registry.Get<std::vector<std::string> >("observables", "Common", "Observables");
    for (unsigned int k = 0; k < observable_ids.size(); ++k)
    {
        std::transform ( observable_ids[k].begin(), observable_ids[k].end(), observable_ids[k].begin(), ::toupper );//the result is stored in observable_ids[k]
        if (t_exp != 0.0)
        {
            std::string sigf(propertycache[17].name);
            std::string mgf(propertycache[18].name);
            std::transform ( sigf.begin() , sigf.end() , sigf.begin() , ::toupper );
            std::transform ( mgf.begin() , mgf.end() , mgf.begin() , ::toupper );
            if ((observable_ids[k] == sigf))
            {
                std::cout<<propertycache[17].name<<" makes no sense for a realtime simulation. Ignoring"<<std::endl;
                continue;
            }
            if ((observable_ids[k] == mgf) && !my_model_is_impurity_model)
            {
                std::cout<<propertycache[18].name<<"is only possible for the impurity models in imaginary-time simulations!"<<std::endl;
                continue;
            }
        }
        bool stop = false;
        unsigned int idx = 0;
        while ((!stop) && (idx < nr_of_known_observables) )
        {
	  const ObservableProperties& props(propertycache[idx]);
            std::string comp(props.name);
            std::transform ( comp.begin() , comp.end() , comp.begin() , ::toupper );
            if (comp == observable_ids[k])
            {
	      if(!props.makes_only_sense_for_impurity || my_model_is_impurity_model)
	      {
		Early_ObservableProperties earlyprop(idx);
		if((props.type == VECTORFUNCTION) || (props.type == FUNCTION))
		{
		  try
		  {
		    if(registry.Get<bool>("observables", props.name, "covariance"))
		    {
		      std::cout<<"Covariance-Matrix requested for "<< props.name<<std::endl;
		      earlyprop.covariance = true;
		    }
		  }
		  catch (Registry_Key_not_found_Exception& e)
		  {
		    std::cout<<e.what()<<std::endl;
		  }
		}
                observablenames.push_back(earlyprop);
                std::cout<<"Adding Observable: "<<propertycache[observablenames.back().idx].name<<std::endl;
                stop = true;
	      }
	      else
		std::cout<<"Not adding "<<propertycache[idx].name<<"! It makes no sense for an impurity model."<<std::endl;
            }
            else
                ++idx;
        }
        if (idx == nr_of_known_observables)
        {
            std::cout<<"Observable not recognized: "<<observable_ids[k]<<std::endl;
        }
    }
}
catch (Registry_Exception& e)
{
    cout<<e.what()<<endl;
    throw;
}
catch (...)
{
    throw;
}

/**
A Helper function to generate the observables using the right data type. The data type depends on the template parameters
*/
template <typename ObsType, typename SignType>
void complex_obs_helper(std::vector<ObservableBase<SignType>*>& observables, const Early_ObservableProperties& earlyprop, const Parameters& curparams)
{
  const ObservableProperties& props(propertycache[earlyprop.idx]);
    switch (props.type)
    {
    case SCALAR:
        observables.push_back(new Observable<ObsType, SignType >(props.name, curparams));
        break;
    case FUNCTION:
        observables.push_back(new Observable<std::valarray<ObsType>, SignType>(props.name, curparams, earlyprop.covariance));
        break;
    case VECTORFUNCTION:
        observables.push_back(new Observable<std::valarray<std::valarray<ObsType> >, SignType>(props.name, props.isk_space ?"k":"r", curparams, earlyprop.covariance));
        break;
    default:
        throw(std::runtime_error("this shouldn't happen... "));
    }
}

//Next follow various helpers that instantiate the complex_obs_helper to distinguish between the various possibilities of getting a complex or real sign
template <class Sign>
void branch_on_sign(bool obscomplex, std::vector<ObservableBase<Sign>*>& observables, const Early_ObservableProperties& earlyprop, const Parameters& curparams);
template <>
void branch_on_sign<double>(bool obscomplex, std::vector<ObservableBase<double>*>& observables, const Early_ObservableProperties& earlyprop , const Parameters& curparams)
{
    if (obscomplex)
    {
        complex_obs_helper<std::complex<double>, double>(observables, earlyprop, curparams);
    }
    else
    {
        complex_obs_helper<double, double>(observables, earlyprop, curparams);
    }
    return;
}

template <>
void branch_on_sign<std::complex<double> >(bool obscomplex, std::vector<ObservableBase<std::complex<double> >*>& observables, const Early_ObservableProperties& earlyprop, const Parameters& curparams)
{
    if (obscomplex)
        complex_obs_helper<std::complex<double>, std::complex<double> >(observables, earlyprop, curparams);
    else
        throw(std::runtime_error("Hey YOU! How should having a complex sign but a real observable work out?"));
    return;
}

template <typename SignType>
void Master::wait_for_all_clients_to_stop(std::list<ClientState<SignType> >& clients, std::ofstream& additionalinfo)
{
    std::cout<<"waiting for clients to stop"<<std::endl;
    bool all_stopped = false;
    while (!all_stopped)
    {
        try
        {
            int num = checkActivity();//see if sth. is happening on any port
	    if (num > 0)
	    {
	      while (datapending())
                {
                    const Client& o = peekNextOriginator();
                    typename std::list<ClientState<SignType> >::iterator it = clients.begin();
                    for (; (it != clients.end()) && it->client != o; ++it);
                    if (it != clients.end())
                    {
		if (peeknextSize(o) == 1)
                {
                    Letter<char> ack(recvfromspecificClient(o));
                    if (ack.msg == STOP_RECEIVED)
                    {
                        std::cout<<"Stop ACK received"<<std::endl;
			it->running = false;
                    }
                }
                else
                {
                    discardmessage(o);//can't be an ack. hence we discard that message
                    std::cout<<"Discarding a message!"<<std::endl;
                }
                    }
                    else
                    {
                        cout<<"A serious error happened! But we don't care!"<<endl;
                        exit(-1);
                    }
                }
	    }
        }
        catch (ConnectionError& e)
        {
            const string msg( e.client.getName() + string(" decapitated."));
            time_t rawtime;
            time ( &rawtime );
            cout<<"["<<ctime(&rawtime)<<"] "<<"Activity on Client Socket but no Data received! Assuming Zombie finally dead!"<<endl;
            cout<<msg<<endl;
            additionalinfo<<msg<<std::endl;
            clients.erase(remove_if(clients.begin(), clients.end(), Client_Eraser<SignType>(e.client)), clients.end());//remove the client from our array
            removeClient(e.client);//remove the client from the net subsystem
        }
        typename std::list<ClientState<SignType> >::iterator it = clients.begin();
        all_stopped = (it->running == false);
        ++it;
        for ( ;it != clients.end(); ++it)
        {
            all_stopped = ((!it->running) && all_stopped);
        }
    }
    std::cout<<"All clients have stopped"<<std::endl;
}

template <typename SignType>
int Master::handleActivity(std::list<ClientState<SignType> >& clients, std::ofstream& additionalinfo)
{
    int num = 0;
    try
    {
        num = checkActivity();//see if sth. is happening on any port
    }
    catch (ConnectionError& e)
    {
        const string msg( e.client.getName() + string(" decapitated."));
	time_t rawtime;
        time ( &rawtime );
        cout<<"["<<ctime(&rawtime)<<"] "<<"Activity on Client Socket but no Data received! Assuming Zombie finally dead!"<<endl;
        cout<<msg<<endl;
        additionalinfo<<msg<<std::endl;
        clients.erase(remove_if(clients.begin(), clients.end(), Client_Eraser<SignType>(e.client)), clients.end());//remove the client from our array
        removeClient(e.client);//remove the client from the net subsystem
        num = 0;
    }
    return num;
}

template <typename SignType>
void Master::run(std::vector<Parameters>::const_iterator curparams, std::list<ClientState<SignType> >& clients)
{
    bool run = true;
    AverageSign<SignType> sign(curparams->binpath);
    std::vector<ObservableBase<SignType>*> observables;//herein we store all observables
    //create the various observables
    for (typename std::vector<Early_ObservableProperties>::const_iterator it = observablenames.begin(); it != observablenames.end(); ++it)
    {
        const bool iscomplexobservable = (curparams->signiscomplex || propertycache[it->idx].iscomplex);
        branch_on_sign<SignType>(iscomplexobservable, observables, *it, *curparams);
    }
    const unsigned int currentsize = observables.front()->bin_nr();
    for (int k = 1; k < static_cast<int>(observables.size()); ++k)
        if (currentsize != observables[k]->bin_nr())
        {
            cout<<"Inconsistent bin sizes detected: bins: "<<observables[k]->bin_nr()<<endl;
            exit(-1);
        }
    ofstream additionalinfo((curparams->path + curparams->idpath + std::string("AdditionalInfo.txt")).c_str(), ios::app);
    time_t rawtime;
    time ( &rawtime );
    additionalinfo<<"beginning to raise undead at: "<<ctime (&rawtime)<<std::endl;
    if ( usetimer || (currentsize < cycles) ||
      (analyzeonly == false)/* && 
      (currentsize < cycles) && !usetimer*/
    )
    {
        cout<<"Waiting for new souls........"<<endl;
        //send the new parameters to all clients
        for (typename std::list<ClientState<SignType> >::iterator client = clients.begin(); client != clients.end(); ++client)
        {
            connectclient(client->client, *curparams);
            time ( &rawtime );
            additionalinfo<<"Taking control of: "<<client->client.getName()<<" at "<<ctime (&rawtime)<<std::endl;
            while (true)
            {
	      std::cout<<"looping1"<<std::endl;
                try
                {
                    checkActivity();//see if sth. is happening on any port
                }
                catch (ConnectionError& e)
                {
                    const string msg( e.client.getName() + string(" decapitated."));
                    time_t rawtime;
                    time ( &rawtime );
                    cout<<"["<<ctime(&rawtime)<<"] "<<"Activity on Client Socket but no Data received! Assuming Zombie finally dead!"<<endl;
                    cout<<msg<<endl;
                    additionalinfo<<msg<<std::endl;
		    if(e.client == client->client)//the disconnected client is the same as the client that we're waiting for
		    {
		      client = clients.erase(client);
		      removeClient(e.client);//remove the client from the net subsystem
		      break;
		    }
		    else
		    {
		      clients.erase(client);
		      removeClient(e.client);//remove the client from the net subsystem
		    }
                }
                if (anyMessagefrom(client->client))
                {
                    if (peeknextSize(client->client) == 1)
                    {
                        Letter<char> ack(recvfromspecificClient(client->client));
                        if (ack.msg == INIT_RECEIVED)
			{
			  std::cout<<"ACK from "<<client->client.getName()<<" received!"<<std::endl;
			  client->running = true;
			  break;
			}
                    }
                    else
                    {
                        discardmessage(client->client);//can't be an ack. hence we discard that message
                    }
                }
            }
        }
        while (run)
        {
            int num = 0;//handleActivity<SignType>(clients, additionalinfo);
            try
            {
                num = checkActivity();//see if sth. is happening on any port
            }
            catch (ConnectionError& e)
            {
                const string msg( e.client.getName() + string(" decapitated."));
                time ( &rawtime );
                cout<<"["<<ctime(&rawtime)<<"] "<<"Activity on Client Socket but no Data received! Assuming Zombie finally dead!"<<endl;
                cout<<msg<<endl;
                additionalinfo<<msg<<std::endl;
                clients.erase(remove_if(clients.begin(), clients.end(), Client_Eraser<SignType>(e.client)), clients.end());//remove the client from our array
                removeClient(e.client);//remove the client from the net subsystem
                continue;
            }
            if (num > 0)//there's some data on the ports!
            {
                if (clientringing())
                {
                    std::cout<<"[Master]: Nr of Clients: "<<clients.size()<<std::endl;
                    const Client& c = addClient();
                    connectclient(c, *curparams);
                    if (wait_for_Init_Ack(c, clients, additionalinfo) == true)
                    {
                        clients.push_back(ClientState<SignType>(c, sign, observables));
                        clients.back().running = true;
                    }
                }
                while (datapending())
                {
                    const Client& o = peekNextOriginator();
                    typename std::list<ClientState<SignType> >::iterator it = clients.begin();
                    for (; (it != clients.end()) && it->client != o; ++it);
                    if (it != clients.end())
                    {
                        it->addnextData(*this);
                    }
                    else
                    {
                        cout<<"A serious error happened! But we don't care!"<<endl;
                        exit(-1);
                    }
                    if (!usetimer)
                    {
                        if (sign.size() >= cycles) run = false;
                    }
                    else
                    {
                        double elapsed = elapsedtime();
                        if (elapsed >= interval) run = false;
                    }
                }
            }
        }
        clearInbox();
        sendtoallclients<unsigned char>(STOP);
	sign.sync();
	for (int k = 0; k < static_cast<int>(observables.size()); ++k)//sync for SuperMUC
	{
	  observables[k]->sync();
	}
        wait_for_all_clients_to_stop(clients, additionalinfo);
    }
    additionalinfo<<"QMC Run finished."<<std::endl;
    //now do error analysis...
    const unsigned int checksize = observables.front()->bin_nr();
    additionalinfo<<"Nr of measurements taken : "<<checksize<<endl;
    for (int k = 1; k < static_cast<int>(observables.size()); ++k)
        if (checksize != observables[k]->bin_nr())
        {
            cout<<"Inconsistent bin sizes detected!"<<endl;
            exit(-1);
        }
#ifdef _OPENMP
    double start = omp_get_wtime();
#endif
    sign.analyzeData(curparams->path + curparams->idpath);
#pragma omp parallel for schedule(dynamic)
    for (int k = 0; k < static_cast<int>(observables.size()); ++k)
        observables[k]->analyzeData(curparams->path + curparams->idpath, sign);
    //tidy up observables
#ifdef _OPENMP
    std::cout<<"Error Analysis took "<< omp_get_wtime() -start<<" seconds."<<std::endl;
#endif
    for (unsigned int k = 0; k < observables.size(); ++k)
        delete observables[k];
    if (curparams->model == IMAG_MODEL_FROM_FILE )
        delete [] curparams->data;
}

Master::~Master()
{
    sendtoallclients<unsigned char>(TERMINATE);
    sleep(2);//the zombies sleep for one second. give them a chance to wake up
    cout<<"successfully terminating"<<endl;
}

int main(int argc, char *argv[])
{
    Master master(argc, argv);
    std::list<ClientState<std::complex<double> > > clientscomplex;
    std::list<ClientState<double> > clients;
    for(std::vector<Parameters>::const_iterator curparams = master.params.begin(); curparams != master.params.end(); ++curparams)
   {
    std::cout<<"starting new simulation"<<std::endl;
    if (curparams->signiscomplex)
        master.run<std::complex<double> >(curparams, clientscomplex);
    else
        master.run<double>(curparams, clients);
    std::cout<<"finished run!"<<std::endl;
    }
    return 0;
}

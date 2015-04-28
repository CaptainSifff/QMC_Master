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
#ifndef MPISERVER_H
#define MPISERVER_H
#define OMPI_SKIP_MPICXX
#include <mpi.h>
#include <vector>
#include <valarray>
#include <complex>
#include <iostream>
#include <stdint.h>
#include "MPITypeDeduction.h"

using namespace std;

class Client
{
public:
    inline string getName() const;
    unsigned int id;
    inline Client(int myid) : id(myid) {}
    inline Client() {}
    inline bool operator==(const Client& a) const throw();
    inline bool operator!=(const Client& a) const throw();
};

#include "MPI_Error.h"

/**
Compares two Clients.
@return true if it is the same socket and the same Client-Name
*/
bool Client::operator==(const Client& a) const throw()
{
    return (a.id == id);
}

bool Client::operator!=(const Client& a) const throw()
{
    return !(*this == a);
}

std::string Client::getName() const
{
    stringstream ss;
    ss<< id;
    return ss.str();
}

#include "NetError.h"

template <typename T>
class Letter
{
public:
    Client sender;
    T msg;
};

class MPIServer
{
public:
    /**
    With this Function you can check whether a client wants to establish a connection
    @return true if a client wants to establish a connection, else false
    */
    inline bool clientringing();
    /**
    This function will return if there is activity on any connection. If not it waits indefinitely(~49 days)
    @return the number of active ports
    */
    inline int checkActivity();
    /**
    Constructor call. It needs argc and argv from the command line, to pass it on to MPI
    @param argc argc as given from the main funtion
    @param argv argv as given from the main function
    */
    inline MPIServer(int argc, char *argv[]);
    inline ~MPIServer();
    /**
    With this you can send sth. to a client.
    Note to self: MPI_Send(void *buf, ... is the MPI standard but as it seems MPI doesn't modify the buffer (which it would be allowed, judging from the signature)
    we can use a const_cast to make it fit...
    */
    template <typename T>
    inline bool send(const Client&, const T&);
    /**
    This function returns a message from any client.
    The content already has the type specified by the template parameter.
    */
    template <typename T>
    inline Letter<T> recvfromanyclient();
    inline Letter<char> recvfromspecificClient(const Client&);
    inline Client peekNextOriginator();
    /**
    This function adds a Client that's waiting for a connection to the internal structures
    @return the added client
    */
    inline const Client& addClient();
    /** This functions returns the size of the next message.
     * @param 
     */
    inline uint peeknextSize(const Client& c);
    /**
    This function removes a client(this event can't happen MPI)
    @return It returns alwys zero right now...
    */
    inline int removeClient(const Client&);
    /**
    This function checks wether there is data waiting on any connection
    @return true if there is data on any connection
    */
    inline bool datapending();
    /**
    With this function you can send some data to all connected clients
    */
    template <typename T>
    inline bool sendtoallclients(T);
    /** This function enables checking for a message from a particular client.
     * @param client The client for which to check.
     * @return true if MPI_Probe finds some message, else false.
     */
    inline bool anyMessagefrom(const Client& client);
    /** This discards the message of the client c
     * @param c the client whose message to discard
     * */
    inline void discardmessage(const Client& c);
    void clearInbox() {} //hopefully MPI takes care of this...
private:
    bool juststarted;
    bool hasdata;
    int MPInrofclients;
    int alreadyaddedclients;
    vector<Client> clients;
    MPI_Status status;///< here we store the state of the MPI_Status variable. It enables fetching things like size.
    bool status_is_uptodate;
    /** checks the status of the internal status variable
    */
    void status_valid();
};

inline void MPIServer::discardmessage(const Client& c)
{
  if(status_is_uptodate)
  {
    int len = peeknextSize(c);
    char* dump = new char[len]; 
    MPI_Recv(dump, len, MPI_CHAR, c.id, status.MPI_TAG, MPI_COMM_WORLD, &status);
    status_is_uptodate = false;
    delete [] dump;
  }
  else
  {
    	  std::cout<<"status not up-to-date in discardmessage !"<<std::endl;
	  exit(-1);
  }
}

uint MPIServer::peeknextSize(const Client& c)
{
    int length;
    if(status_is_uptodate)
        MPI_Get_count(&status, MPI_CHAR, &length);//hopefully this works for all datatypes...
	else
	{
	  std::cout<<"status not up-to-date in peeknextsize !"<<std::endl;
	  exit(-1);
	}
    return length;
}

bool MPIServer::anyMessagefrom(const Client& client)
{
  int retval;
  MPI_Iprobe(client.id, MPI_ANY_TAG, MPI_COMM_WORLD, &retval, &status);
  if (!retval) status_is_uptodate = false;
    else
        status_is_uptodate = true;
  return retval;
}

void MPIServer::status_valid()
{
    if (!status_is_uptodate) throw MPI_Error("[MPIServer] MPI_Status not up to date!!");
}

template <typename T>
bool MPIServer::sendtoallclients(T var)
{
    bool success = true;
    for (vector<Client>::const_iterator it = clients.begin(); it != clients.end(); ++it)
    {
        success = (success && send(*it, var));
    }
    return success;
}

int MPIServer::checkActivity()
{
    if (juststarted) return true;//the user code has to connect all the clients first
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);//this call blocks until a matching msg is pending
    hasdata = true;
    status_is_uptodate = true;
    return true;
}

Client MPIServer::peekNextOriginator()
{
    return Client(status.MPI_SOURCE);
}

bool MPIServer::datapending()
{
    int flag;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
    if (!flag) status_is_uptodate = false;
    else
        status_is_uptodate = true;
    return flag;
}

bool MPIServer::clientringing()
{
    return juststarted;
}

template <typename T>
inline bool MPIServer::send(const Client& client, const T& var)
{
    return MPI_Send(const_cast<T*>(&var), 1, mpi::MPIType<T>::get(), client.id, 314, MPI_COMM_WORLD) != 0;
}

template <>
inline bool MPIServer::send<std::string>(const Client& client, const std::string& var)
{
    return MPI_Send(const_cast<char*>(var.c_str()), var.size(), MPI_CHAR, client.id, 314, MPI_COMM_WORLD) != 0;
}

template <>
inline bool MPIServer::send<std::valarray<char> >(const Client& client, const std::valarray<char>& var)
{
    return MPI_Send(const_cast<char*>(&(var[0])), var.size(), MPI_CHAR, client.id, 314, MPI_COMM_WORLD) != 0;
}

template <>
inline Letter<valarray<double> > MPIServer::recvfromanyclient<valarray<double> >()
{
    status_valid();
    Letter<valarray<double> > retval;
    int length;
    MPI_Get_count(&status, MPI_DOUBLE, &length);
    retval.msg.resize(length);
    MPI_Recv(&(retval.msg[0]), length, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
    retval.sender.id = status.MPI_SOURCE;
    status_is_uptodate = false;
    return retval;
}

template <>
inline Letter<valarray<complex<double> > > MPIServer::recvfromanyclient<valarray<complex<double> > >()
{
//this implementation uses the OpenMPI MPI_DOUBLE_COMPLEX define
    status_valid();
    Letter<valarray<complex<double> > > retval;
    int length;
    int err = MPI_Get_count(&status, MPI_DOUBLE_COMPLEX, &length);
    if (err != 0)
    {
        cout<<"Error while receiving the length of the array! Aborting."<<endl;
        exit(-1);
    }
    retval.msg.resize(length);
    MPI_Recv(&(retval.msg[0]), length, MPI_DOUBLE_COMPLEX, /*status.MPI_SOURCE*/MPI_ANY_SOURCE, MPI_ANY_TAG/*status.MPI_TAG*/, MPI_COMM_WORLD, &status);
    retval.sender.id = status.MPI_SOURCE;
    status_is_uptodate = false;
    return retval;
}

template <>
inline Letter<std::valarray<valarray<complex<double> > > > MPIServer::recvfromanyclient<std::valarray<valarray<complex<double> > > >()
{
//this implementation uses the OpenMPI MPI_DOUBLE_COMPLEX define
    status_valid();
    typedef valarray<complex<double> > Function;
    Letter<std::valarray<Function> > retval;
    uint32_t funnr;
    MPI_Recv(&funnr, 1, MPI_UNSIGNED, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);//first get the number of function points. This is the message that checkactivity reacted upon.
    retval.msg.resize(funnr);
    //now get the actual function points. we need to poll again for matching data
    MPI_Probe(status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);//ask MPI politely to provide the number of the transmitted function points
    int length;
    if (MPI_Get_count(&status, MPI_DOUBLE_COMPLEX, &length) != 0)
    {
        std::cout<<"Error while receiving the length of the array! Aborting."<<std::endl;
        exit(-1);
    }
    if (length == MPI_UNDEFINED)
    {
        std::cout<<"the length of the buffer was not defined!"<<std::endl;
        exit(-1);
    }
    std::complex<double>* temp = new std::complex<double>[length];
    MPI_Recv(temp, length, MPI_DOUBLE_COMPLEX, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
    unsigned int functionpoints = length/funnr;
    for (unsigned int k = 0; k < funnr; ++k)
        retval.msg[k].resize(functionpoints);
    for (unsigned int k = 0; k < funnr; ++k)
        for (unsigned int j = 0; j < functionpoints; ++j)
            retval.msg[k][j] = temp[k*functionpoints + j];
    delete [] temp;
    retval.sender.id = status.MPI_SOURCE;
    status_is_uptodate = false;
    return retval;
}

template <>
inline Letter<std::valarray<valarray<double> > > MPIServer::recvfromanyclient<std::valarray<valarray<double> > >()
{
    status_valid();
    typedef valarray<double> Function;
    Letter<std::valarray<Function> > retval;
    uint32_t funnr;
    MPI_Recv(&funnr, 1, MPI_UNSIGNED, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);//first get the number of function points. This is the message that checkactivity reacted upon.
    retval.msg.resize(funnr);
    //now get the actual function points. we need to poll again for matching data
    MPI_Probe(status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);//ask MPI politely to provide the number of the transmitted function points
    int length;
    if (MPI_Get_count(&status, MPI_DOUBLE, &length) != 0)
    {
        std::cout<<"[MPIServer] Error while receiving the length of the array! Aborting."<<std::endl;
        exit(-1);
    }
    if (length == MPI_UNDEFINED)
    {
        std::cout<<"[MPIServer] The length of the buffer was not defined!"<<std::endl;
        exit(-1);
    }
    double* temp = new double[length];
    MPI_Recv(temp, length, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
    unsigned int functionpoints = length/funnr;
    for (unsigned int k = 0; k < funnr; ++k)
        retval.msg[k].resize(functionpoints);
    for (unsigned int k = 0; k < funnr; ++k)
        for (unsigned int j = 0; j < functionpoints; ++j)
            retval.msg[k][j] = temp[k*functionpoints + j];
    delete [] temp;
    retval.sender.id = status.MPI_SOURCE;
    status_is_uptodate = false;
    return retval;
}

template <typename T>
inline Letter<T> MPIServer::recvfromanyclient()
{
    status_valid();
    Letter<T> retval;
    MPI_Recv(&(retval.msg), 1, mpi::MPIType<T>::get(), status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
    retval.sender.id = status.MPI_SOURCE;
    status_is_uptodate = false;
    return retval;
}

inline Letter<char> MPIServer::recvfromspecificClient(const Client& c)
{
    status_valid();
    Letter<char> retval;
    int length;
    MPI_Get_count(&status, MPI_CHAR, &length);
    MPI_Recv(&(retval.msg), length, MPI_CHAR, c.id, status.MPI_TAG, MPI_COMM_WORLD, &status);
    retval.sender.id = c.id;
    status_is_uptodate = false;
    return retval;
}

const Client& MPIServer::addClient()
{
    clients.push_back(Client(1 + alreadyaddedclients) );
    alreadyaddedclients++;
    if (alreadyaddedclients >= MPInrofclients)
    {
        juststarted = false;
    }
    return clients.back();
}

int MPIServer::removeClient(const Client&)
{
    return 0;//to suppress gcc warning
}

MPIServer::MPIServer(int argc, char *argv[]): juststarted(true), hasdata(false), alreadyaddedclients(0), status_is_uptodate(false)
{
    int numprocs, myid;
    int flag = 0;
    MPI_Initialized(&flag);
    if (!flag)
    {
        if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
        {
            throw MPI_Error("[MPIServer] Error initializing MPI!");//fake exception... MPI shall per default abort if an error happens...
        }
    }
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPInrofclients = numprocs - 1;
    return;
}

MPIServer::~MPIServer()
{
    MPI_Finalize();
}
#endif

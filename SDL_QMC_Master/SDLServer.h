/***************************************************************************
 *   Copyright (C) 2009 - 2011 by Florian Goth   *
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
#include <sstream>
#include <valarray>
#include <stdexcept>
#include "NetCommonBase.h"
#include "utility.h"
#include "NetServer.h"
#include "NetError.h"
#include "charconverters.h"

using namespace std;

template <typename T>
class Letter
{
public:
    Client sender;
    T msg;
    typedef T value_type;
};

class SDLServer : public SDL_Net_Server::NetServer
{
public:
//    bool clientringing();//supplied by the NetServer
    inline int checkActivity();
    /**
    Constructor call. At the moment we don't use these parameters, as they have no meaning for the SDL
    @param argc argc as given from the main funtion
    @param argv argv as given from the main function
    */
    inline SDLServer(int& argc, char *argv[]);
    inline virtual ~SDLServer();
    template <typename T>
    inline bool send(const Client&, const T&);
    /**
    This function returns a message from any client.
    The content already has the type specified by the template parameter.
    @return A Letter class that contains
    */
    template <typename T>
    inline Letter<T> recvfromanyclient();
    inline Letter<char> recvfromspecificClient(const Client&);
//    const Client& addClient();//supplied by the Netserver
//    int removeClient(const Client&);//supplied by the Netserver
//    bool datapending();//supplied by the NetServer
    /**
     * A function to send some message to all clients that are attached
     */
    template <typename T>
    inline bool sendtoallclients(T);
    inline void discardmessage(const Client& c);
private:
};

void SDLServer::discardmessage(const Client& c)
{
  getNextLetter(c);
}

inline Letter<char> SDLServer::recvfromspecificClient(const Client& c)
{
  Letter<char> retval;
  SDL_Net_Server::Letter letter = getNextLetter(c);
  retval.sender = letter.getsender();
  retval.msg = *(letter.getmsg());
  return retval;
}

template <class T>
struct Priv_SDL
{
    typedef T value_type;
    static inline void recvpack(SDL_Net_Server::Letter& letter, value_type& retval)
    {
        CH<T> temp(letter.getmsg());
        retval = temp;
    }
    static inline void sendpack(const T& arg, Message& m)
    {
        CH<T> cd(arg);
        m.alloc(cd.size());
        for (unsigned int j = 0; j < cd.size(); ++j)
            m.msg[j] = cd.accessbytes(j);
    }
};

template <class T>
struct Priv_SDL<std::valarray<T> >
{
    typedef std::valarray<T> value_type;
    static inline void recvpack(SDL_Net_Server::Letter& letter, value_type& retval)
    {
        const unsigned int len = letter.getlen();
        const unsigned int elems = len/sizeof(T);
        retval.resize(elems);
        for (unsigned int k = 0; k < elems; ++k)
        {
            CH<T> temp(&(letter.getmsg()[k*sizeof(T)]));
            retval[k] = temp;
        }
    }
    static inline void sendpack(const value_type& arg, Message& m)
    {
        m.alloc(arg.size()*sizeof(CH<T>));
        for (unsigned int k = 0; k < arg.size(); ++k)
        {
            CH<T> cd(arg[k]);
            for (unsigned int j = 0; j < cd.size(); ++j)
                m.msg[k * sizeof(CH<T>) + j] = cd.accessbytes(j);
        }
    }
};

template <class T>
struct Priv_SDL<std::valarray<std::valarray<T> > >
{
    typedef std::valarray<std::valarray<T> > value_type;
    static inline void recvpack(SDL_Net_Server::Letter& letter, value_type& retval)
    {
        const uint32_t nroffunctions(CH<uint32_t>(letter.getmsg()));//the first four bytes hold the number of functions
        const char* const msgbegin = &(letter.getmsg()[sizeof(uint32_t)]);//get the new offset
        retval.resize(nroffunctions);
        const unsigned int len = letter.getlen() - sizeof(uint32_t);//get the new index
        const unsigned int functionpoints = len/sizeof(T)/nroffunctions;
        for (unsigned int k = 0; k < nroffunctions; ++k)
            retval[k].resize(functionpoints);
        for (unsigned int j = 0; j < nroffunctions; ++j)
            for (unsigned int k = 0; k < functionpoints; ++k)
            {
                CH<T> temp(&(msgbegin[(j * functionpoints + k) * sizeof(T)]));
                retval[j][k] = temp;
            }
    }
    static inline void sendpack(const value_type& arg, Message& m)
    {
        const uint32_t nrf = arg.size();
        CH<uint32_t> ui(nrf);
        unsigned int len = arg[0].size();
        m.alloc(sizeof(uint32_t) + nrf * len * sizeof(CH<T>));//the first four bytes hold the number of functions
        for (unsigned int j = 0; j < sizeof(uint32_t); ++j)
            m.msg[j] = ui.accessbytes()[j];
        for (unsigned int i = 0; i < nrf; ++i)
        {
            for (unsigned int k = 0; k < len; ++k)
            {
                CH<T> cd(arg[i][k]);
                for (unsigned int j = 0; j < cd.size(); ++j)
                    m.msg[i *len * sizeof(CH<T>)  + k * sizeof(CH<T>) + j + sizeof(uint32_t)] = cd.accessbytes(j);
            }
        }
    }
};

template <class T>
inline Letter<T> SDLServer::recvfromanyclient()
{
    Letter<T> retval;
    SDL_Net_Server::Letter letter = getNextLetter();
    retval.sender = letter.getsender();
    Priv_SDL<T>::recvpack(letter, retval.msg);
    return retval;
}

template <>
inline bool SDLServer::sendtoallclients<unsigned char>(unsigned char arg)
{
    char a = arg;
    SDL_Net_Server::NetServer::sendtoallClients(Message(&a, 1));
    return true;
}

template <>
inline bool SDLServer::send<char>(const Client& c, const char& arg)
{
    char t = arg;
    SDL_Net_Server::NetServer::send(c, Message(&t, 1));
    return true;
}

template <typename T>
inline bool SDLServer::send(const Client& c, const T& arg)
{
    Message m;
    Priv_SDL<T>::sendpack(arg, m);
    SDL_Net_Server::NetServer::send(c, m);
    m.free();
    return true;
}

int SDLServer::checkActivity()
{
    return SDL_Net_Server::NetServer::checkActivity(-1);//we wait indefinitely, to emulate the mpi like behaviour
}

SDLServer::SDLServer(int& argc, char *argv[])
{
    if (SDL_Init(0) == -1)//Initialize SDL
    {
        const string ret("Error in initializing the SDL library!");
        SDL_Quit();
        throw SDL_Error(string("[SDLServer] ") + ret);
    }
    if (SDLNet_Init() == -1)
    {
        const string ret("Error in initializing the SDL_net library!");
        SDLNet_Quit();
        SDL_Quit();
        throw SDL_Error(string("[SDLServer] ") + ret);
    }
    std::string lastarg(argv[argc-1]);
    std::string prelastarg(argv[argc-2]);
    if(prelastarg == "-p")
    {//some dirty hack...
      std::istringstream temp(lastarg);
      uint16_t port;
      if (!(temp >> port) )
	throw invalid_argument("parsing of port failed!");
//    = static_cast<uint_16t>(atoi(lastarg.c_str()));
      argc -= 2;
      this->init(port);
    }
    else
    this->init(31415);
    return;
}

SDLServer::~SDLServer()
{
    SDLNet_Quit();
    SDL_Quit();
    return;
}

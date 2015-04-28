/***************************************************************************
 *   Copyright (C) 2008,2009 by Florian Goth   *
 *   CaptainSifff@gmx.de   *
 *                                                                         *
 *   Permission is hereby granted, free of charge, to any person obtaining *
 *   a copy of this software and associated documentation files (the       *
 *   "Software"), to deal in the Software without restriction, including   *
 *   without limitation the rights to use, copy, modify, merge, publish,   *
 *   distribute, sublicense, and/or sell copies of the Software, and to    *
 *   permit persons to whom the Software is furnished to do so, subject to *
 *   the following conditions:                                             *
 *                                                                         *
 *   The above copyright notice and this permission notice shall be        *
 *   included in all copies or substantial portions of the Software.       *
 *                                                                         *
 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,       *
 *   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF    *
 *   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. *
 *   IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR     *
 *   OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, *
 *   ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR *
 *   OTHER DEALINGS IN THE SOFTWARE.                                       *
 ***************************************************************************/
#ifndef NETSERVER_H
#define NETSERVER_H
#include "NetCommonBase.h"
#include "Client.h"
#include "NetError.h"
#include <algorithm>
#include <vector>
#include <queue>
#include <fstream>
#include <cstring>

using namespace std;

namespace SDL_Net_Server
{
class Letter
{
public:
    inline const Client& getsender() const throw();
    inline char* const& getmsg() const throw();
    inline Letter(Client&, const Message&);
    inline Letter(const Letter&);
#ifdef HAS_RVALUE_REFERENCES
    inline Letter(Letter&& rhs) : sender(rhs.sender), payload(rhs.payload), len(rhs.len)
    {
      rhs.payload = NULL;
    }
    inline Letter& operator=(Letter&& rhs)
    {
      if(this != &rhs)
      {
	sender = rhs.sender;
	len = rhs.len;
	delete [] payload
	payload = rhs.payload;
	rhs.payload = NULL;
      }
      return *this;
    }
#endif
    inline Letter& operator=(const Letter& rhs)
    {
      if(this != &rhs)
      {
	sender = rhs.sender;
	len = rhs.len;
	delete [] payload;
	payload = new char[len];
	memcpy(payload, rhs.payload, len);
      }
      return *this;
    }
    inline ~Letter();
    inline unsigned int getlen() const throw();
private:
    /*const */Client sender;
    char* payload;
    unsigned int len;
};

Letter::Letter(const Letter& rhs) : sender(rhs.sender), payload(new char[rhs.len]), len(rhs.len)
{
memcpy(payload, rhs.payload, len);
return;
}

Letter::~Letter()
{
delete [] payload;
}

const Client& Letter::getsender() const throw()
{
    return sender;
}

unsigned int Letter::getlen() const throw()
{
    return len;
}

char* const& Letter::getmsg() const throw()
{
    return payload;
}

Letter::Letter(Client& c, const Message& m) : sender(c), payload(new char[m.len]), len(m.len)
{
    memcpy(payload, m.msg, m.len);
}

struct ClientInfo
{
    Client client;
    uint32_t ping;
    uint32_t pong;
    bool pingpong;
    ClientInfo() : ping(0), pong(0), pingpong(false) {}
    ClientInfo& operator=(const ClientInfo&) = default;
    ClientInfo(Client a) : client(a), ping(0), pong(0), pingpong(false) {}
};

class NetServer : protected NetCommonBase<TCPsocket>
{
public:
    inline int sendtoallClients(const Message);
    inline unsigned int send(TCPsocket, const Message&);
    inline unsigned int send(const Client&, const Message&);
    inline int init(uint16_t);
    inline int checkActivity(uint32_t) throw(ConnectionError);
    inline bool clientringing() throw();
    inline const Client& addClient();
    inline int removeClient(const Client&);
    inline NetServer();
    inline virtual ~NetServer();
    inline bool datapending() throw();
    inline Letter getNextLetter();
    inline Letter getNextLetter(const Client&);
    inline const string& findClient(TCPsocket);
    inline const Client& findClient(const string&);
    inline vector<string> listClientNames();
    inline float getPing(const Client&) throw();
    inline Client peekNextOriginator();
    inline uint peeknextSize() const {return inbox.front().getlen();}
    inline uint peeknextSize(const Client& c) const;
    inline char *const peekNextContent() const {return inbox.front().getmsg();}
    inline bool anyMessagefrom(const Client&);
    typedef deque<Letter> InBoxType;
    void clearInbox(){inbox.clear();}
protected:
    inline unsigned int nrofClients() const;
private:
    InBoxType inbox;
    IPaddress myIP;
    TCPsocket mysocket;
    string myname;
    vector<ClientInfo> Clients;
    SDLNet_SocketSet socketset;//the socketset where we watch for activity
    inline int setupSocketset();
//    deque<Letter> inbox;
    bool clientisringing;
    inline void processPong(ClientInfo&) throw();
    inline void pingall() throw();
    inline void ping(Client&) throw();
    inline ClientInfo& findbyClient(const Client&) throw(NetError);
    inline InBoxType::const_iterator findnextClientmessage(const Client& c) const {
      InBoxType::const_iterator it = inbox.begin();
       for(; (it != inbox.end()) && (it->getsender() != c); ++it )
	 ;
	 return it;
    }
    inline InBoxType::iterator findnextClientmessage(const Client& c) {
      InBoxType::iterator it = inbox.begin();
       for(; (it != inbox.end()) && (it->getsender() != c); ++it )
	 ;
	 return it;
    }
};

uint NetServer::peeknextSize(const Client& c) const
{
    InBoxType::const_iterator it = findnextClientmessage(c);//this might segfault if there is noc actual message...
    return it->getlen();
}

bool NetServer::anyMessagefrom(const Client& c)
{
  InBoxType::iterator it;
  for(it = inbox.begin(); (it != inbox.end()) && (it->getsender() != c); ++it )
    ;
  return inbox.end() != it;
}

Client NetServer::peekNextOriginator()
{
  return inbox.front().getsender();
}

void NetServer::ping(Client& c) throw()
{
    char ping[] = {'S', 'E', 'R', 'V', 'E', 'R', '_', 'P', 'I', 'N', 'G', 0};
    send(c, Message(ping, strlen(ping)));
}

float NetServer::getPing(const Client& cl) throw()
{
    try
    {
        ClientInfo& ci = findbyClient(cl);
        cout<<ci.pong<<" "<<ci.ping<<endl;
        return static_cast<float>(ci.pong - ci.ping)/1000.0f;
    }
    catch (NetError& e)
    {
        cout<<e.what()<<endl;
    }
    return 0;
}

void NetServer::pingall() throw()
{
    for (vector<ClientInfo>::iterator it = Clients.begin() ; it != Clients.end(); ++it)
    {
        it->ping = SDL_GetTicks();
        it->pingpong = true;
        ping(it->client);
    }
    return;
}

vector<string> NetServer::listClientNames()
{
    vector<string> ret;
    for (vector<ClientInfo>::const_iterator it = Clients.begin(); it != Clients.end(); ++it)
        ret.push_back(it->client.getName());
    return ret;
}

unsigned int NetServer::nrofClients() const
{
    return Clients.size();
}

const string& NetServer::findClient(TCPsocket sock)
{
    vector<ClientInfo>::const_iterator it = Clients.begin();
    for (; (it != Clients.end()) && (it->client.getSocket() != sock); ++it) {};
    return it->client.getName();
}

const Client& NetServer::findClient(const string& str)
{
    vector<ClientInfo>::const_iterator it = Clients.begin();
    for (; (it != Clients.end()) && (it->client.getName() != str); ++it);
    return it->client;
}

NetServer::~NetServer()
{
    //clean up the connections to the remaining clients
    while(!Clients.empty())
      removeClient(Clients.back().client);
    if(socketset)
      SDLNet_FreeSocketSet(socketset);
    if(mysocket)
      SDLNet_TCP_Close(mysocket);
    return;
}

int NetServer::sendtoallClients(const Message m)
{
    int ret = 0;
    for (vector<ClientInfo>::const_iterator it = Clients.begin(); it != Clients.end(); ++it)
    {
        ret = send(it->client.getSocket(), m);
        if (ret != static_cast<int>(m.len))
            ret = -1;
    }
    return ret;
}

unsigned int NetServer::send(const Client& client, const Message& m)
{
    return NetCommonBase<TCPsocket>::send(client.getSocket(), m);
}

unsigned int NetServer::send(TCPsocket s , const Message& m)
{
    return NetCommonBase<TCPsocket>::send(s, m);
}

bool NetServer::datapending() throw()
{
    return  !inbox.empty();
}

NetServer::NetServer() : clientisringing(false)
{
    return;
}

int NetServer::setupSocketset()
{
    socketset = SDLNet_AllocSocketSet(Clients.size() +1 );
    if (!socketset)
    {
        cout<<"SDLNet_AllocSocketSet: "<<SDLNet_GetError()<<endl;
        return -1;
    }
    //Attach Server
    SDLNet_TCP_AddSocket(socketset, mysocket);
    //Attach Clients
    for ( vector<ClientInfo>::iterator it = Clients.begin(); it != Clients.end(); ++it)
        SDLNet_TCP_AddSocket(socketset, it->client.getSocket());
    return 0;
}

/**
Adds the Client that's waiting for a Connection
@return A reference to the freshly connected client
*/
const Client& NetServer::addClient()
{
    TCPsocket sock = SDLNet_TCP_Accept(mysocket);
    if (sock != NULL)
    {
        IPaddress* clientaddress = SDLNet_TCP_GetPeerAddress(sock);
        const char* hostname = SDLNet_ResolveIP(clientaddress);
        string hoststring( hostname == NULL ? "N/A": hostname);
        cout<<"Connection from Client "<<hoststring<<"("<<
        convertIPtoString(clientaddress->host)<<":"<<clientaddress->port<<") accepted"<<endl;
        Message username;
        if ((username = recv(sock,username)).msg == NULL)
        {
            cout<<"Error in receiving username!"<<endl;
            SDLNet_TCP_Close(sock);
            throw UserNameError("Error in receiving username!");
        }
        string un(username.msg);
        username.free();
        Clients.push_back(ClientInfo(Client(sock, un)));
//Update Socket Sets
        if (socketset)
            SDLNet_FreeSocketSet(socketset);
        if (setupSocketset() != 0)
            throw NetError("Unknown Network Error");
    }
    else
    {
        SDLNet_TCP_Close(sock);
        throw ErrorConnecting("[Net] [TCP] Couldn't accept Connection");
    }
    clientisringing = false;
    pingall();//Update ping Times
    return Clients.back().client;
}

class compareClients : public std::unary_function<Client, bool>
{
public:
    inline bool operator()(const ClientInfo& arg) const 
    {
        return arg.client == comp;
    }
    compareClients(const Client& a) : comp(a) {}
private:
    const Client& comp;
};

ClientInfo& NetServer::findbyClient(const Client& c) throw(NetError)
{
    vector<ClientInfo>::iterator it = find_if(Clients.begin(), Clients.end(), compareClients(c));
    if (it == Clients.end())
    {
        cout<<"Client "<<c.getName()<<" not found!"<<endl;
        throw(NetError("Client " + c.getName() + " not found!"));
    }
    return *it;
}

int NetServer::removeClient(const Client& arg)
{
    vector<ClientInfo>::iterator it = find_if(Clients.begin(), Clients.end(), compareClients(arg));
    if (it == Clients.end())
    {
        cout<<"Client "<<arg.getName()<<" not found!"<<endl;
        return -1;
    }
    SDLNet_TCP_Close(it->client.getSocket());
    Clients.erase(it);
    SDLNet_FreeSocketSet(socketset);
    return setupSocketset();
}

bool NetServer::clientringing() throw()
{
    return clientisringing;
}

Letter NetServer::getNextLetter(const Client& c)
{
    deque<Letter>::iterator it = findnextClientmessage(c);
    Letter retval(*it);
    inbox.erase(it);
    return retval;
}

Letter NetServer::getNextLetter()
{
    Letter retval(inbox.front());
    inbox.pop_front();
    return retval;
}

void NetServer::processPong(ClientInfo& cl ) throw()
{
    cl.pingpong = false;
    cl.pong = SDL_GetTicks();
    return;
}

int NetServer::checkActivity(uint32_t t) throw(ConnectionError)
{
    uint32_t t1 = SDL_GetTicks();
    int numactivity = SDLNet_CheckSockets(socketset, t);
    if (numactivity == -1)
    {
        printf("SDLNet_CheckSockets: %s\n",SDLNet_GetError());
        return -1;
    }
    int events = numactivity;
    if (SDLNet_SocketReady(mysocket))
    {
        clientisringing = true;
        events--;
    }
    for (vector<ClientInfo>::iterator it = Clients.begin(); (it != Clients.end()) && (events > 0); it++)
        if (SDLNet_SocketReady(it->client.getSocket()))
        {
            Message temp;
            try
            {
                temp = recv(it->client.getSocket(), temp);
            }
            catch (TimeOutError& e)
            {
                vector<ClientInfo>::const_iterator it2 = Clients.begin();
                for (; (it2 != Clients.end()) && (it2->client.getSocket() != e.sock); ++it2) {};
                throw ConnectionError(string(e.what()),it2->client);
            }
            catch (NetError& e)
            {
                cout<<e.what()<<endl;
                throw;
            }
            if ((temp.msg != NULL) && (temp.len > 0))
            {
                static char const pong[] = {'C','L', 'I', 'E', 'N','T','_','P','O','N','G', 0};
                if ((temp.len == strlen(pong)) && (strncmp(temp.msg, pong, strlen(pong)) == 0))
                {
                    numactivity--;
                    processPong(*it);
                    cout<<"Ping Time: "<<getPing(it->client)<<endl;
                }
                else
                    inbox.push_back(Letter(it->client, temp));
            }
            else
                throw( ConnectionError("[Net] Client seems to be dead...", it->client));
            events--;
            temp.free();
        }
    uint32_t t2 = SDL_GetTicks();
    unsigned int dt = t2-t1;
    if ((numactivity == 0) && (dt < t))//We processed a Pong Event, but we shouldn't return yet
        return checkActivity(t - dt);
    return numactivity;
}

int NetServer::init(uint16_t port)
{
    //Set up for using a server socket
    if (SDLNet_ResolveHost(&myIP,NULL,port) == -1)
    {
        const string ret("SDLNet_ResolveHost: " + string(SDLNet_GetError()));
        throw SDL_Error(string("[SDLServer] ") + ret);
    }
    // output the IP address nicely
    cout<<"IP Address: "<<convertIPtoString(myIP.host)<<endl;

    // resolve the hostname for the IPaddress
    const char* host = SDLNet_ResolveIP(&myIP);

    // print out the hostname we got
    printf("Hostname   : %s\n",(host ? host :"N/A"));
    // open the server socket
    mysocket = SDLNet_TCP_Open(&myIP);
    if (!mysocket)
    {
        const string ret("SDLNet_TCP_Open: " + string(SDLNet_GetError()));
        throw SDL_Error(string("[SDLServer] ") + ret);
    }
    //Add the server to the socketset
    socketset = SDLNet_AllocSocketSet(1);
    if (!socketset)
    {
        const string ret("SDLNet_AllocSocketSet: " + string(SDLNet_GetError()));
        throw SDL_Error(string("[SDLServer] ") + ret);
    }
    SDLNet_TCP_AddSocket(socketset, mysocket);
    return 0;
}
}
#endif

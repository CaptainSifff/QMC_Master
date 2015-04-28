/***************************************************************************
 *   Copyright (C) 2008 by Florian Goth   *
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
#ifndef NETERROR_H
#define NETERROR_H
#include <string>
#include <exception>

class NetError : public std::exception
{
public:
    NetError(std::string a): msg(a) {}
    ~NetError() throw() {}
    const char* what() const throw()
    {
        return msg.c_str();
    }
private:
    std::string msg;
};

class SubSystemError : public NetError
{
public:
    SubSystemError(std::string a) : NetError(a) {}
private:
};

class ErrorConnecting : public NetError
{
public:
    ErrorConnecting(std::string a) : NetError(a) {}
private:
};

class UserNameError : public NetError
{
public:
    UserNameError(std::string a) : NetError(a) {}
private:
};

class ConnectionError : public NetError
{
public:
    ConnectionError(std::string a, Client c) : NetError(a), client(c){}
    ~ConnectionError() throw() {}
    Client client;
private:
};
#endif

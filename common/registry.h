/***************************************************************************
 *   Copyright (C) 2005-2009 by Florian Goth   *
 *   Captainsifff@gmx.de   *
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
#ifndef REGISTRY_H
#define REGISTRY_H

/**
@author Florian Goth
*/
#include <algorithm>
#include <map>
#include <string>
#include <vector>
#include <stdexcept>
#include <sstream>

/**
The basic Exception thrown by the registry. It is derived from the STL logic_error exception
*/
class Registry_Exception : public std::logic_error
{
public:
    explicit Registry_Exception(const std::string& what_arg) throw(): std::logic_error(what_arg) {}
private:
};

/**
The Exception when the desired key is not found
*/
class Registry_Key_not_found_Exception : public Registry_Exception
{
public:
    explicit Registry_Key_not_found_Exception(const std::string& what_arg) throw(): Registry_Exception(what_arg) {}
private:
};

class Registry_Block_Data_not_found_Exception : public Registry_Key_not_found_Exception
{
public:
    const std::string block_key;
    explicit Registry_Block_Data_not_found_Exception(const std::string& key_value) throw(): Registry_Key_not_found_Exception(std::string("Block Key not found: ") + key_value), block_key(key_value) {}
    ~Registry_Block_Data_not_found_Exception() throw() {}
private:
};

class Registry_Block_not_found_Exception : public Registry_Key_not_found_Exception
{
public:
    const std::string block;
    explicit Registry_Block_not_found_Exception(const std::string& key_value) throw(): Registry_Key_not_found_Exception(std::string("Block not found: ") + key_value), block(key_value) {}
    ~Registry_Block_not_found_Exception() throw() {}
private:
};

class Registry_cfgfile_not_found_Exception : public Registry_Key_not_found_Exception
{
public:
    const std::string cfgfile;
    explicit Registry_cfgfile_not_found_Exception(const std::string& key_value) throw(): Registry_Key_not_found_Exception(std::string("cfgfile not found: ") + key_value), cfgfile(key_value) {}
    ~Registry_cfgfile_not_found_Exception() throw() {}
private:
};

/**
this class contains the contents of a [BLOCK] structure in a config file.
*/
class Block
{
    //In Config Files those [BLOCK] thingies
    std::map<std::string , std::string > Block_Data;
    std::string BlockName;
    std::vector<std::string> Keys;
    std::map<std::string , std::string >::size_type NrOfKeys;
    void push_back_Key(std::string);
public:
    Block() {}
    Block(std::string blockname , const std::vector<std::string>& block);
    inline unsigned int size(void) const;
    const std::vector<std::string> GetKeys() const
    {
        return Keys;
    }
    /**
    get access to the contents of a block
    */
    std::string& operator[](const std::string& Bn)
    {
        return Block_Data[Bn];
    }
    std::string& find(const std::string& key) throw(Registry_Block_Data_not_found_Exception)
    {
        std::map<std::string , std::string >::iterator it = Block_Data.find(key);
        if (it == Block_Data.end()) throw(Registry_Block_Data_not_found_Exception(key));
        return it->second;
    }
    const std::string& find(const std::string& key) const throw(Registry_Block_Data_not_found_Exception)
    {
        std::map<std::string , std::string >::const_iterator it = Block_Data.find(key);
        if (it == Block_Data.end()) throw(Registry_Block_Data_not_found_Exception(key));
        return it->second;
    }
    ~Block()
    {
    }
};

/**
this class contains the contents of a single config file
*/
class cfile
{
    std::string filename;
    std::map<std::string , Block> cfgfile;
public:
    cfile()
    {}
    Block& operator[](std::string a )
    {
        return cfgfile[a];
    }
    Block& find(const std::string& key) throw(Registry_Block_not_found_Exception)
    {
        std::map<std::string , Block>::iterator it = cfgfile.find(key);
        if ( it == cfgfile.end()) throw(Registry_Block_not_found_Exception(key));
        return it->second;
    }
    const Block& find(const std::string& key) const throw(Registry_Block_not_found_Exception)
    {
        std::map<std::string , Block>::const_iterator it = cfgfile.find(key);
        if ( it == cfgfile.end()) throw(Registry_Block_not_found_Exception(key));
        return it->second;
    }
    /**
    constructor that initializes this with the contents of a Config-File
    */
    cfile(std::string& file);
};

/**
this holds together all the contents of the Configuration Directory and provides access via the Get() template
*/
class RegistryDB
{
    std::map<std::string , cfile> Reg;
private:
public:
    /**
    this initializes the registry
    @param arg the Directory that contains all the files the registry should contain
    */
    int Init(const std::string& cfgDir);
    /**
    this is the constructor for the registry
    @param arg the Directory that contains all the files the registry should contain
    */
    RegistryDB(const std::string& arg);
    RegistryDB()
    {}
    ~RegistryDB()
    {}
    Block GetBlock(std::string file , std::string bloc)
    {
        return Reg[file][bloc];
    }
    /**
    function to get a value from the registry. The template parameter determines to which type to convert the key.
    @param file the file in which to look
    @param block under which block is the value
    @param key for which key to look
    @return the requested key
    */
    template < typename T >
    inline T Get(const std::string& file, const std::string& block, const std::string& key) const throw(Registry_Key_not_found_Exception);
    template <typename T>
    inline T set(const std::string& file, const std::string& block, const std::string& key, T val);
};

/**
the basic template for doing the conversion between strings and the requested type
We use the C++ stringstreams thus we benefit from all overloads that are already provided by C++
*/
template < typename A >
struct GetTrait
{
    static inline A Convert(std::string arg)
    {
        A ret;
        std::stringstream sstream(arg);
        sstream >> ret;
        return ret;
    }
};

template <>
struct GetTrait<bool>
{
    static bool Convert(std::string arg)
    {
        std::transform ( arg.begin() , arg.end() , arg.begin() , ::toupper );
        if ( arg == "TRUE")
            return true;
        return false;
    }
};

template < typename A >
struct SetTrait
{
    static std::string convert(A arg)
    {
        std::stringstream ss;
        ss << arg << '\0';
        std::string os(ss.str());
        return os;
    }
};

template < typename T >
inline T RegistryDB::Get(const std::string& file, const std::string& block, const std::string& Key) const throw(Registry_Key_not_found_Exception)
{
    typedef GetTrait<T> Trait;
    std::map<std::string , cfile>::const_iterator tempkey = Reg.find(file);
    if (tempkey == Reg.end()) throw(Registry_cfgfile_not_found_Exception(std::string("File Key not found :") + file));
    std::string temp(tempkey->second.find(block).find(Key));
    return Trait::Convert(temp);
}

template <typename T>
inline T RegistryDB::set(const std::string& file, const std::string& block, const std::string& key, T val)
{
    typedef SetTrait<T> Trait;
    Reg[file][block][key] = Trait::convert(val);
    return val;
}

inline unsigned int Block::size(void) const
{
    return Block_Data.size();
}
#endif

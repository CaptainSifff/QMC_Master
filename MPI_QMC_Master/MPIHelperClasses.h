#ifndef MPI_HELPER_CLASSES_H
#define MPI_HELPER_CLASSES_H

#include <stdint.h>
#include <mpi.h>
#include <valarray>
#include <vector>
#include <complex>
namespace mpi
{
	/**
	 * \brief MPI type deduction class
	 * Specializations of this class should be used to derive the MPI Datatype
	 * from a template parameter.
	 */
	template <class T> class MPIType;

	template <> struct MPIType<double> {static MPI_Datatype get(void){return MPI_DOUBLE;} };
	template <> struct MPIType<float> {static MPI_Datatype get(void){return MPI_FLOAT;} };
        template <> struct MPIType<std::complex<double> > {static MPI_Datatype get(void){return MPI_DOUBLE_COMPLEX;} };
        template <> struct MPIType<std::complex<float> > {static MPI_Datatype get(void){return MPI_COMPLEX;} };
	template <> struct MPIType<uint32_t> {static MPI_Datatype get(void){return MPI_UNSIGNED;} }; // FIXME: Unsure about size of MPI_UNSIGNED
	template <> struct MPIType<int32_t> {static MPI_Datatype get(void){return MPI_INT;} }; // FIXME: Unsure about size of MPI_INT
	template <> struct MPIType<char> {static MPI_Datatype get(void){return MPI_CHAR;} };
	template <> struct MPIType<unsigned char> {static MPI_Datatype get(void){return MPI_UNSIGNED_CHAR;} };
	template <> struct MPIType<long unsigned int> {static MPI_Datatype get(void){return MPI_UNSIGNED_LONG;} };
}


//class MPIServer;

template <class C, class MPINode>
struct BcastHelper
{
	static void bcast(MPINode * const , C & var, int root)
	{
		MPI_Bcast(&var,1,mpi::MPIType<C>::get(), root, MPI_COMM_WORLD);
	}
};



template <class C, class Server_, class Client_>
struct SendHelper//<C, Server_, Client_>
{//generic case
		static bool send(const Client_& client, Server_*, const C& var, const int & tag)
		{
			return MPI_Send(const_cast<C*>(&var), 1, mpi::MPIType<C>::get(), client.id, tag, MPI_COMM_WORLD) != 0; 
		}
};

// Valarray sending as array
template <class C, class Client_ , class Server_>
struct SendHelper<std::valarray<C>, Server_, Client_>
{
		static bool send(const Client_& client, Server_ *const, const std::valarray<C>& var, const int & tag)
		{
			return MPI_Send(const_cast<C*>(&var[0]), var.size(), mpi::MPIType<C>::get(), client.id, tag, MPI_COMM_WORLD) != 0; 
		}
};

template <class C, class Client_>
C RecvHelper<C, Client_>::recv(Client_*, const int & tag)
{
	MPI_Status status;
	C ret;
	MPI_Recv(&ret, 1, mpi::MPIType<C>::get(), MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
	return ret;
}

// Valarray receiving as array
template <class C, class Client_>
struct RecvHelper<std::valarray<C>, Client_>
{
		static std::valarray<C> recv(Client_*, const int & tag)
		{
			MPI_Status status;
			int length;
			MPI_Probe( MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status); // Blocking probe!
			MPI_Get_count(&status, mpi::MPIType<C>::get(), &length);
			std::valarray<C> ret(length);
			MPI_Recv(&ret[0], length, mpi::MPIType<C>::get(), MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
			return ret;
		}
};

// Vector sending as array
template <class C, class Client, class Server_>
struct SendHelper<std::vector<C>, Server_, Client>
{
		static bool send(const Client& client, Server_*, const std::vector<C>& var, const int & tag)
		{
			return MPI_Send(const_cast<C*>(&var[0]), static_cast<int>(var.size()), mpi::MPIType<C>::get(), client.id, tag, MPI_COMM_WORLD) != 0; 
		}
};

// Vector receiving as array
template <class C, class Client_>
struct RecvHelper<std::vector<C>, Client_>
{
		static std::vector<C> recv(Client_*, const int & tag)
		{
			MPI_Status status;
			int length;
			MPI_Probe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status); // Blocking probe!
			MPI_Get_count(&status, mpi::MPIType<C>::get(), &length);
			std::vector<C> ret(length);
			MPI_Recv(&ret[0], length, mpi::MPIType<C>::get(), MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
			return ret;
		}
};

//
//
// Sending from client to server
//
//
template<class C, class Client>
struct SendToServerHelper
{//generic case
	static void sendtoserver(Client * const client, const C & var, const int & tag)
	{
		MPI_Send( const_cast<C*>(&var) , 1, mpi::MPIType<C>::get(), 0, tag, MPI_COMM_WORLD);
	}
};

template<class C, class Client>
struct SendToServerHelper<std::valarray<C>, Client>
{//valarray sending
	static void sendtoserver(Client * const client, const std::valarray<C> & var, const int & tag)
	{
		MPI_Send( const_cast<C*>(&var[0]), static_cast<int>(var.size()), mpi::MPIType<C>::get(), 0, tag, MPI_COMM_WORLD);
	}
};

template<class C, class Client>
struct SendToServerHelper<std::vector<C>, Client>
{//vector sending
	static void sendtoserver(Client * const client, const std::vector<C> & var, const int & tag)
	{
		MPI_Send( const_cast<C*>(&var[0]), static_cast<int>(var.size()), mpi::MPIType<C>::get(), 0, tag, MPI_COMM_WORLD);
	}
};


//
//
// Receiving messages from client
//
//
// FIXME: This can be definitely done by only _one_ kind of RecvHelpers
// independent of the direction Server->Client or Client->Server.
template<typename T> class Letter;

template<class C, class Server, class Client>
struct RecvFromClientHelper
{
	static Letter<C> recv(const Client & client, Server * const server, const int & tag)
	{
		Letter<C> ret;
		MPI_Status status;
		MPI_Recv(&ret.msg,1,mpi::MPIType<C>::get(), client.id, tag, MPI_COMM_WORLD, &status);  
		ret.sender.id=status.MPI_SOURCE;
		return ret;
	}
};

template<class Scalar, class Server, class Client>
struct RecvFromClientHelper<std::valarray<Scalar>, Server, Client>
{// valarray receiving
	static Letter<std::valarray<Scalar> > recv(const Client & client, Server * const server, const int & tag)
	{
		Letter<std::valarray<Scalar> > retval;
		MPI_Status status;
		int length;
		MPI_Probe(client.id, tag, MPI_COMM_WORLD, &status);

		MPI_Get_count(&status, mpi::MPIType<Scalar>::get(), &length);
		retval.msg.resize(length);
		MPI_Recv(&(retval.msg[0]), length, mpi::MPIType<Scalar>::get() , status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
		retval.sender.id = status.MPI_SOURCE;
		return retval;

	}
};


template<class Scalar, class Server, class Client>
struct RecvFromClientHelper<std::vector<Scalar>, Server, Client>
{// vector receiving
	static Letter<std::vector<Scalar> > recv(const Client & client, Server * const server, const int & tag)
	{
		Letter<std::vector<Scalar> > retval;
		MPI_Status status;
		int length;
		MPI_Probe(client.id, tag, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, mpi::MPIType<Scalar>::get(), &length);
		retval.msg.resize(length);
		MPI_Recv(&(retval.msg[0]), length, mpi::MPIType<Scalar>::get() , status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
		retval.sender.id = status.MPI_SOURCE;
		return retval;

	}
};




#endif

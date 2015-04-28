#include <complex>
namespace mpi
{
/**
 * \brief MPI type deduction class
 * Specializations of this class should be used to derive the MPI Datatype
 * from a template parameter.
 */
template <class T> struct MPIType;

template <> struct MPIType<double>
{
    static MPI_Datatype get(void) {return MPI_DOUBLE;}
};
template <> struct MPIType<float>
{
    static MPI_Datatype get(void) {return MPI_FLOAT;}
};
template <> struct MPIType<std::complex<double> >
{
    static MPI_Datatype get(void) {return MPI_DOUBLE_COMPLEX;}
};
template <> struct MPIType<std::complex<float> >
{
    static MPI_Datatype get(void) {return MPI_COMPLEX;}
};
template <> struct MPIType<uint32_t>
{
    static MPI_Datatype get(void) {return MPI_UNSIGNED;}
}; // FIXME: Unsure about size of MPI_UNSIGNED
template <> struct MPIType<int32_t>
{
    static MPI_Datatype get(void) {return MPI_INT;}
}; // FIXME: Unsure about size of MPI_INT
template <> struct MPIType<char>
{
    static MPI_Datatype get(void) {return MPI_CHAR;}
};
template <> struct MPIType<unsigned char>
{
    static MPI_Datatype get(void) {return MPI_UNSIGNED_CHAR;}
};
template <> struct MPIType<long unsigned int>
{
    static MPI_Datatype get(void) {return MPI_UNSIGNED_LONG;}
};
}

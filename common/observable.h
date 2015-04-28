#ifndef OBSERVABLE_H
#define OBSERVABLE_H
#include "filevector.h"
#include "AverageSign.h"
#include "Master.h"

template<class Source>
inline std::string toString(const Source& Value)
{
    std::stringstream ss;
    ss << Value;
    std::string os(ss.str());
    return os;
}

template <class SignType>
class ObservableBase
{
public:
    inline ObservableBase(const std::string s) : myname(s) {}
    virtual void commit() = 0;
    virtual void cacheData(Master&) = 0;
    virtual void analyzeData(const std::string&, AverageSign<SignType>&) = 0;
    virtual unsigned int bin_nr() = 0;
    virtual void sync() = 0;
    virtual ~ObservableBase() {}
    const std::string myname;
};

template <class ObsT, class SignType>
class Observable : public ObservableBase<SignType>
{
public:
    template <class T, class S> friend class Error_Analysis;
    typedef ObsT ObsType;
    void commit();
    void cacheData(Master&);
    void analyzeData(const std::string&, AverageSign<SignType>&);
    Observable(const std::string&, const Parameters&);
    void sync() {fv.sync();}
    unsigned int bin_nr()
    {
        return fv.size();
    }
private:
    Observable(const Observable&);
    Observable& operator=(Observable&);
    ObsType cache;
    FileVector<ObsType> fv;
};

template <class ObsT, class SignType>
Observable<ObsT, SignType>::Observable(const std::string& m, const Parameters& p) : ObservableBase<SignType>(m)
{
    fv.load(p.binpath + m);//FIXME: we load from binpath + name of observable
}

template <class ObsT, class SignType>
void Observable<ObsT, SignType>::commit()
{
    fv.push_back(cache);
}

template <class ObsT, class SignType>
void Observable<ObsT, SignType>::cacheData(Master& m)
{
    cache = m.recvfromanyclient<ObsType>().msg;
}

//provide a specialization for valarrays
template <class ObsT, class SignType>
class Observable<std::valarray<ObsT>, SignType> : public ObservableBase<SignType>
{
public:
    template <class T, class S> friend class Error_Analysis;
    typedef std::valarray<ObsT> ObsType;
    void commit();
    void cacheData(Master&);
    void analyzeData(const std::string&, AverageSign<SignType>&);
    Observable(const std::string&, const Parameters&, bool);
    unsigned int bin_nr()
    {
        return fv.size();
    }
    void sync() {fv.sync();}
private:
    Observable(const Observable&);
    Observable& operator=(Observable&);
    const uint32_t functionpoints;
    ObsType cache;
    FileVector<ObsType> fv;
    const double delta_s;
    bool covariance;
};

template <class ObsT, class SignType>
Observable<std::valarray<ObsT>, SignType>::Observable(const std::string& m, const Parameters& p, bool cv = false) : ObservableBase<SignType>(m), functionpoints(
  m.compare(0, 10, "KondoCloud") == 0 ? p.Nb*2*p.Nx:
  p.functionpoints
), cache(functionpoints), delta_s(p.delta_s), covariance(cv)
{
    fv.load(p.binpath + m, functionpoints);//FIXME: we load from binpath + name of observable
}

template <class ObsT, class SignType>
void Observable<std::valarray<ObsT>, SignType>::commit()
{
    fv.push_back(cache);
}

//provide a specialization for valarrays<valarrays>, aka tensor functions
template <class ObsT, class SignType>
class Observable<std::valarray<std::valarray<ObsT> >, SignType> : public ObservableBase<SignType>
{
public:
    template <class T, class S> friend class Error_Analysis;
    typedef std::valarray<ObsT> Function;
    typedef std::valarray<Function> ObsType;
    void commit();
    void cacheData(Master&);
    void analyzeData(const std::string&, AverageSign<SignType>&);
    Observable(const std::string&, const std::string&, const Parameters&, bool);
    unsigned int bin_nr();
    void sync() {fv.sync();}
private:
    Observable(const Observable&);
    Observable& operator=(Observable&);
    std::string indexname;
    const unsigned int tensorindices;
    ObsType cache;
    FileVector<ObsType> fv;
    const uint32_t& functionpoints;
    const double delta_s;
    std::string name;
    bool covariance;
};

template <class ObsT, class SignType>
void Observable<std::valarray<std::valarray<ObsT> >, SignType>::cacheData(Master& m)
{
    cache = (m.recvfromanyclient<std::valarray<Function> >()).msg;
    return;
}

template <class ObsT, class SignType>
unsigned int Observable<std::valarray<std::valarray<ObsT> >, SignType>::bin_nr()
{
    return fv.size();
}

template <class ObsT, class SignType>
void Observable<std::valarray<std::valarray<ObsT> >, SignType>::commit()
{
  fv.push_back(cache);
}

template <class ObsT, class SignType>
Observable<std::valarray<std::valarray<ObsT> >, SignType>::Observable(const std::string& m, const std::string& in, const Parameters& p, bool cv = false) : ObservableBase<SignType>(m), indexname(in),
        tensorindices(
        ( p.is_Impurity_model?
        (p.model == SIAM?p.N
        :
        (((p.model == KONDO_IMP_TI) && (m.substr(0, 24) == "LocalBathGreensfunctions"))? p.Nx*p.Nb*2: p.N)
        )
        :p.N)
        ),
        cache(tensorindices),
        functionpoints(p.functionpoints),
        delta_s(p.delta_s), name(m), covariance(cv)
{
    fv.load(p.binpath + m, functionpoints, tensorindices);//FIXME: we load from binpath + name of observable
    for (unsigned int k = 0; k < tensorindices; ++k)
    {
        cache[k].resize(functionpoints);
    }
    return;
}
#endif

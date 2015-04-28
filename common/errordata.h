# ifndef MONTE_CARLO_TOOLS_ERRORDATA_H // include guard
# define MONTE_CARLO_TOOLS_ERRORDATA_H
#include <valarray>
#include <iostream>
namespace mc_analysis
{

template <typename R> class errordata
{
private:
    R mean;
    R error;
    R bias;///< the bias of the jackknife procedure. This is e.g useful for estimating the quality of the binsize. WARNING: not all functions of the mc_analysis class set this
public:
    /**
    A constructor for the error data of a scalar observable
    @param m this sets the estimation for the mean value of observable
    @param e this set the estimation for the standard error of the observable
    @param b this set the estimate for the bias
    */
    errordata(const R& m, const R& e, const R& b) : mean(m), error(e), bias(b) {}
    /**
    Two argument version of the constructor. The bias is set to zero.
    @param m this sets the estimation for the mean value of observable
    @param e this set the estimation for the standard error of the observable
    */
    errordata(const R& m, const R& e) : mean(m), error(e), bias(0) {}
    void set_mean(R mean__)
    {
        mean=mean__;
    };
    void set_error(R error__)
    {
        error=error__;
    };
    R get_mean(void) const
    {
        return mean;
    }
    R get_error(void) const
    {
        return error;
    }
    /**
    An accessor to the bias of the observable
    @return the bias of the observable
    */
    R get_bias() const
    {
        return bias;
    }
    typedef R value_type;
//    size_t size() const {return }
};

template <typename R>
inline std::ostream& operator<<(std::ostream& out, const errordata<R>& arg )
{
    out<<arg.get_mean()<<" +- "<<arg.get_error()<<std::endl;
    out<<"bias: "<<arg.get_bias();
    return out;
}

// Specialization for valarray:
template <typename Scalar> class errordata<std::valarray<Scalar> >
{
public:
    typedef std::valarray<Scalar> CovType;
    void setCov(std::valarray<Scalar>& rhs) {
      cov.resize(rhs.size());
      cov = rhs;      
    }
    const CovType& getCov() const {return cov;}
    /**
    A constructor for the error data of a scalar observable
    @param m this sets the estimation for the mean value of observable
    @param e this set the estimation for the standard error of the observable
    @param b this set the estimate for the bias
    */
    errordata(const std::valarray<Scalar>& m, const std::valarray<Scalar>& e, const std::valarray<Scalar>& b) : mean(m), error(e), bias(b) {}
    void set_mean(std::valarray<Scalar> mean_ )
    {
        mean.resize( mean_.size() );
        mean=mean_;
    }
    void set_error(std::valarray<Scalar> error_ )
    {
        error.resize( error_.size() );
        error=error_;
    }
    std::valarray<Scalar> get_mean(void) const
    {
        return mean;
    }
    std::valarray<Scalar> get_error(void) const
    {
        return error;
    }
    /**
    An accessor to the bias of the observable
    @return the bias of the observable
    */
    std::valarray<Scalar> get_bias() const
    {
        return bias;
    }
    errordata<std::valarray<Scalar> > & operator=(const errordata<std::valarray<Scalar> > & errd__)
    {
        this->mean.resize(errd__.mean.size() );
        this->mean=errd__.mean;
        this->error.resize(errd__.error.size() );
        this->error=errd__.error;
	this->cov.resize(errd__.cov.size());
	this->cov=errd__.cov;
        return *this;
    }
    private:
    std::valarray<Scalar> mean;
    std::valarray<Scalar> error;
    std::valarray<Scalar> bias;
    CovType cov;
};

template <typename R> class errordata_signed
{
private:
    R mean;
    R error;
    R sign;
public:
    inline void set_mean(R mean__)
    {
        mean=mean__;
    };
    inline void set_error(R error__)
    {
        error=error__;
    };
    inline void set_sign(R sign__)
    {
        sign=sign__;
    };
    inline R get_mean(void)
    {
        return mean;
    };
    inline R get_error(void)
    {
        return error;
    };
    inline R get_sign(void)
    {
        return sign;
    };
};


}

# endif // include guard

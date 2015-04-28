# ifndef MONTE_CARLO_TOOLS_SCALAR_OBSERVABLE_H  // include guard
# define MONTE_CARLO_TOOLS_SCALAR_OBSERVABLE_H

#include <vector>
#include <valarray>
#include <cmath>

# include "errordata.h"
# include "analysis_common.h"

namespace mc_analysis
{
	class ExampleFunctor
	{
		public:
			typedef	double res_t;
			res_t operator()(double input__){return exp(input__);};
	};

	template <typename S, class Cont = std::vector<S> > class Scalar_Observable : public Cont
		/**
		 * 	\brief This class represents a scalar observable and can be used for the calculation of the mean and
		 * 	error of Monte Carlo data.
		 */
	{
		private:
                        /**
                        As the Scalar_Observable in itself is now a Container all access to the original data points has now to go through Cont::operator[](), and Cont::size()
                        */
			std::valarray<S> data; ///< this valarray stores the data points after the rebinning was applied
			unsigned int num_bins;
			unsigned int binsize;
                        using Cont::size;//make the size function of the base class inaccesible to the outside world
                        using Cont::operator[];
		public:
			inline Scalar_Observable(void);
                        inline virtual ~Scalar_Observable(void);
			inline unsigned int get_num_data_points(void){return Cont::size();};
			inline S mean(void) const;
			inline S sd(void) const;
			inline S standard_error(void) const;
			inline errordata<S> jackknife(void);

			inline void rebin_data(unsigned int binsize__ = 1);
			inline void fractional_rebinning(unsigned int num_bins__);
			inline void print_data(void);
			inline void clear(void);

			template <typename Func> 
                        inline errordata<typename Func::res_t> jack_eval(Func& f);
			inline void operator<<(S bin);
                        using Cont::value_type;
                        /**
                        these access operators are handy for the jackknife analysis and provide access to the binned data
                        @param k the index of the Element
                        @return a reference to the element
                        */
                        inline S& operator[] (std::size_t k) {return data[k];}
                        inline const S& operator[] (std::size_t k) const {return data[k];}
                        inline std::size_t size() const {return data.size();}
	};

	/* Implementation of the class member functions */
       	template <typename S, class Cont> 
        Scalar_Observable<S, Cont>::~Scalar_Observable(void)
	{
            mc_analysis_common::CallClear<Cont>::destroy(static_cast<Cont&>(*this));
        }

	template <typename S, class Cont> 
        Scalar_Observable<S, Cont>::Scalar_Observable(void) : num_bins(0), binsize(0)
	{}
        /**
        This function strips the observable class of all its content. It is intended for use in memory-sparse situations.
        */
	template <typename S, class Cont> 
        void Scalar_Observable<S, Cont>::clear(void)
	{
            data.resize(1);
            num_bins = 0;
            binsize = 0;
            return;
	}

	template <typename S, class Cont> 
        S Scalar_Observable<S, Cont>::mean(void) const
	/**
	 *	\brief This function calculates the arithmetic average of a sample
	 */
	{
		if (num_bins > 0)
		{
			S sum(data[0]);//initialize with first element
			for(unsigned int i = 1; i < num_bins; i++) sum+=data[i];	
			return sum/static_cast<S>(num_bins);
		}
		else
		{
			std::cerr << "Scalar_Observable::mean(): error: no data present! \n";
			return 0;
		}
	}

	template <typename S, class Cont>
        S Scalar_Observable<S, Cont>::sd(void) const
	/**
	 *	\brief This function calculates the standard deviation of the sample
	 */
	{
		S sum = 0.;
		S mean=this->mean();
		for(unsigned int i = 0; i < num_bins; i++) sum += (data[i] - mean)*(data[i] - mean);
		return std::sqrt( sum/static_cast<S> (num_bins - 1.0) );
	}

	template<typename S, class Cont>
        S Scalar_Observable<S, Cont>::standard_error(void) const
	{
		return sd()/std::sqrt( static_cast<S>(num_bins) );
	}

	template<typename S, class Cont> 
        void Scalar_Observable<S, Cont>::rebin_data(unsigned int binsize__)
		/**
		 *	\brief This function has to be called before using the data analysis functions. Rebinning is mandatory
		 *	as it fills the data vector. The original data remains intact in the data_orig vector
		 *	\param binsize__ Specifies, how many data points should be averaged to give a new bin
		 */
	{
                if(1 == binsize__ )
                {//no rebinning -> just copy over
                    data.resize(Cont::size());
                    for(unsigned int k = 0; k < Cont::size(); ++k)
                        data[k] = Cont::operator[](k);
                    this->binsize = 1;
                    this->num_bins = Cont::size();
                }
                else
                {
		   this->binsize = binsize__;
		   this->num_bins = Cont::size()/binsize__;
		   data.resize(num_bins);
		   for(unsigned int i = 0; i < this->num_bins; i++)
		   {
			S sum(0.);
			for(unsigned int j = 0; j < this->binsize; j++)
			{
				sum += Cont::operator[](i*this->binsize + j);
			}
			data[i] = sum/static_cast<S>(this->binsize);
		   }
                }
                return;
	}

	template<typename S, class Cont>
        void Scalar_Observable<S, Cont>::fractional_rebinning(unsigned int num_bins__)
	{
		S prefactor;
		prefactor= static_cast<S>(num_bins__)/static_cast<S>(Cont::size());
		data.resize(num_bins__);	
		this->binsize=0; //here, the binsize would be a fraction, therefore, we can't use the unsigned int
		this->num_bins=num_bins__;
		int n_bins=num_bins__;

		int rest=0;
		int counter=0; // running index
		int total_weight=0;
		for(unsigned int i=0; i<num_bins__; i++)
		{
			data[i]=0.;
			total_weight=0;
			data[i]+= (n_bins-rest)/static_cast<double>(n_bins) * Cont::operator[](counter);
			counter++;
			total_weight+=(n_bins-rest);
			while(total_weight < Cont::size())
			{
				if( (total_weight + n_bins) > Cont::size() )
				{
					rest = Cont::size() - total_weight;
					data[i]+= rest/static_cast<double>(n_bins) * Cont::operator[](counter);
					total_weight = Cont::size();
				}
				else
				{
					data[i] += Cont::operator[](counter);
					total_weight+=n_bins;
					counter++;
				}
			}
			data[i]*=prefactor;
		}
	}

	template<typename S, class Cont>
        void Scalar_Observable<S, Cont>::print_data(void)
	{
		for(unsigned int i=0;i<data.size();i++)
		{
			std::cout<< data[i] << "\t";
		}
		std::cout<<std::endl;
	}

	template <typename S, class Cont>
        void Scalar_Observable<S, Cont>::operator<<(S bin)
	{
		this->push_back(bin);
	}

	template <typename S, class Cont> 
        errordata<S> Scalar_Observable<S, Cont>::jackknife(void)
		/**
		 *	\brief this function performs a simple jackknife analysis for the data
		 *	it calculates the jackknife mean and the standard error.
		 */
	{
		// For this special case, the Jackknife mean is the same as the arithmetic average of the data
		errordata<S> errd;
		errd.set_mean(this->mean());
		// The Jackknife error for this case also reduces to the standard error:
		errd.set_error(this->standard_error());
		return errd;		
	}


	template <typename S, class Cont> 
        template <typename Func> 
        errordata<typename Func::res_t> Scalar_Observable<S, Cont>::jack_eval(Func & f)
		/**
		 * \brief Calculate jackknife error of the observable
		 * This function calculates the mean and error of the data after applying a function
 		 * specified in a Function object.
		 * \param F is a functor defined by the class Functor. You need to overload the operator() as well as
		 * specify the return type of the functor Func::res_t. 
		 */
	{
		// Generate Jackknife bins
		//typename Func::res_t jackdata[this->num_bins];//C99-Style dynamic arrays
		typename Func::res_t* jackdata = new typename Func::res_t[this->num_bins];
		S mean = this->mean();

		for(unsigned int i = 0; i < this->num_bins; i++)
		{
			S x_J;
			x_J = 1./static_cast<S>(this->num_bins - 1) *(this->num_bins *mean - data[i]);
			jackdata[i] = f(x_J);	
		}

		typename Func::res_t jackmean(0.);
		for(unsigned int i=0; i< this->num_bins;i++)
		{
			jackmean += jackdata[i];
		}
		jackmean *= 1.0/static_cast<typename Func::res_t>(this->num_bins);

		typename Func::res_t fun_mean(0.), fun_error(0.);
		fun_mean= f(mean);
		for(unsigned int i=0; i< this->num_bins; i++)
		{
			fun_error += (jackmean - jackdata[i])*(jackmean - jackdata[i]);
		}
		fun_error *= static_cast<typename Func::res_t>(this->num_bins -1)/static_cast<typename Func::res_t>(this->num_bins);
		fun_error = std::sqrt(fun_error);
                delete [] jackdata;

		errordata<typename Func::res_t> errd;
		errd.set_mean(fun_mean);
		errd.set_error(fun_error);
		return errd;
	}
}

#undef GCC_VERSION
# endif  // include guard

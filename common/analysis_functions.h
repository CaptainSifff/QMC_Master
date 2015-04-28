# ifndef MC_ANALYSIS_ANALYSIS_FUNCTIONS_H
# define MC_ANALYSIS_ANALYSIS_FUNCTIONS_H // include guard

template <class T>
struct ToScalar
{
  typedef T RetType;
};

template <class T>
struct ToScalar<std::valarray<T> >
{
  typedef typename ToScalar<T>::RetType RetType;
};

namespace mc_analysis
{
        /**
        This function calculates the expectation value of the samples stored in the container cont
        @param cont the container that stores the samples
        @return the average value of the stored samples
        */
        template<class T>
        typename T::value_type mean(const T& cont)
        {
                 typename T::value_type retval(cont[0]);
                 for (unsigned int k = 1; k < cont.size(); ++k)
                        retval += cont[k];
                 return retval/static_cast<typename ToScalar<T>::RetType >(cont.size());
        }

	template <typename Scalar> Scalar mean_signed(const std::vector<Scalar> & signdata__, const std::vector<Scalar> & sign__)
	{
		Scalar average_sign(0.), average_data(0.);
		average_data=mean(signdata__);
		average_sign=mean(sign__);
		return average_data/average_sign;		
	}

	template <typename Scalar> errordata_signed<Scalar> jackknife_signed(const std::vector<Scalar> & signdata__, const std::vector<Scalar> & sign__)
		/**
		 *	\brief Calculates the jackknife error, mean and average sign of the data.
		 *	\param signdata__ a vector containing the measured data already multiplied by the sign in 
		 *			  the Monte Carlo simulation
		 *	\param sign__	  a vector containing the signs. Note that the sizes of the two vectors have
		 *			  to match!
		 */
	{
		errordata_signed<Scalar> errd;
		Scalar average_sign, average_data;
		average_data=mean(signdata__);
		average_sign=mean(sign__);
		Scalar length= static_cast<Scalar>(signdata__.size());
		std::vector<Scalar> jackdata(signdata__.size());
		// Fill Jackknife bins
		for(unsigned int i=0; i< jackdata.size(); i++)
		{
			jackdata[i] = (length*average_data - signdata__[i])/(length*average_sign - signdata__[i]);	
		}
		Scalar jackmean(0.);
		for(unsigned int i=0; i< signdata__.size();i++)
		{
			jackmean += jackdata[i];
		}
		jackmean *= 1.0/length;

		Scalar mean, error, sign;
		mean= average_data/average_sign;
		for(unsigned int i=0; i< signdata__.size(); i++)
		{
			error+= (jackmean - jackdata[i])*(jackmean - jackdata[i]);
		}
		error *= (length -1.)/length;
		error = std::sqrt(error);

		errd.set_mean(mean);
		errd.set_error(error);
		errd.set_sign(average_sign);
		return errd;
	}


	template <typename Scalar, typename Func> errordata_signed<typename Func::res_t> jack_eval_signed(const std::vector<Scalar> & signdata__, const std::vector<Scalar> & sign__, const Func & f) 
	{
		std::vector<typename Func::res_t> jackdata(signdata__.size());
		errordata_signed<typename Func::res_t> errd;
		Scalar average_sign, average_data;
		average_data=mean(signdata__);
		average_sign=mean(sign__);
		Scalar length= static_cast<Scalar>(signdata__.size());

		for(unsigned int i=0; i<signdata__.size(); i++)
		{
			Scalar x_J; 
			x_J=(length*average_data - signdata__[i])/(length*average_sign - signdata__[i]);
			jackdata[i] = f(x_J);	
		}

		typename Func::res_t jackmean(0.);
		for(unsigned int i=0; i< signdata__.size();i++)
		{
			jackmean += jackdata[i];
		}
		jackmean *= 1.0/static_cast<typename Func::res_t>(signdata__.size());

		typename Func::res_t fun_mean(0.), fun_error(0.);
		fun_mean= f(average_data/average_sign);
		for(unsigned int i=0; i< signdata__.size(); i++)
		{
			fun_error+= (jackmean - jackdata[i])*(jackmean - jackdata[i]);
		}
		fun_error *= static_cast<typename Func::res_t>(signdata__.size() -1)/static_cast<typename Func::res_t>(signdata__.size());
		fun_error = std::sqrt(fun_error);

		errd.set_mean(fun_mean);
		errd.set_error(fun_error);
		errd.set_sign(average_sign);
		return errd;
	
	}


	template <typename Scalar> Scalar autocorrelation_function(const std::vector<Scalar> data, const unsigned int displacement)
		/**
		 * 	\brief Calculate the normalized Autocorrelation function A(displacement). It is defined by:
		 *
		 * 	A(disp)= ( <data[i] data[i+disp]> - <data[i]>^2 )/( <data[i]^2> - <data[i]>^2 )
		 * 	\param data	vector containing the data. The size of the vector should be much larger 
		 * 			than the displacement.
		 * 	\param displacement
		 */
	{
		Scalar sum_disp(0.), sum_mean(0.), sum_sq(0.);	
		for(unsigned int i=0; i< data.size()-displacement; i++)
		{
			sum_disp += data[i]*data[i+displacement];
			sum_mean += data[i];
			sum_sq += data[i]*data[i];	
		}	
		Scalar length=static_cast<Scalar>(data.size() - displacement);
		return (sum_disp - 1.0/(length)*sum_mean*sum_mean )/( sum_sq - 1.0/(length)* sum_mean*sum_mean);
	}

	

}


# endif // include guard

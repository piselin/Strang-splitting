// scalar types
typedef double real_t;
typedef std::complex<real_t> complex_t;
typedef int dim_t;

template<int R, int C>
using CMatrix = Eigen::Matrix<complex_t,R,C>;

template<int R, int C>
using RMatrix = Eigen::Matrix<real_t,R,C>;
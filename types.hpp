// Disclaimer:
// These types are take from libWaveBlocks
// https://github.com/WaveBlocks/libwaveblocks
// to ensure compatibility

// scalar types
using real_t = double;
using complex_t = std::complex<real_t>;
using dim_t = int;

template<int R, int C>
using CMatrix = Eigen::Matrix<complex_t,R,C>;

template<int R, int C>
using RMatrix = Eigen::Matrix<real_t,R,C>;
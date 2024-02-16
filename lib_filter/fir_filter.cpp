//#include <inttypes.h>
//#include <string.h>
//#include <vector>
//
//#include "../../FifoBuffer.h"
//#include "../iir_filter.h"
//#include "iir_filter.h"
//
//class fir_filter: public i_fir_filter
//{
//    fifo_buffer<double> x_ring_;
//    fifo_buffer<double> y_ring_;
//    double n[5];
//    double d[5];
//    size_t nr_coefficients_;
//
//public:
//    fir_filter(const double *n, const double *d, size_t nr_coefficients);
//    double filter(double x);
//    double filter_opt(double x);
//    void init_history_values(double x, int nr_samples = 2000);
//    virtual ~fir_filter() = default;
//};
//
//fir_filter::fir_filter(const double *an, const double *ad, size_t nr_coefficients)
//    : x_ring_(0),
//      y_ring_(0),
//      nr_coefficients_(nr_coefficients)
//{
//    memcpy(n, an, nr_coefficients * sizeof(double));
//    memcpy(d, ad, nr_coefficients * sizeof(double));
//}
//
//double fir_filter::filter(double x)
//{
//    if (x_ring_.size() >= (nr_coefficients_ - 1))
//        return filter_opt(x);
//    else
//    {
//        x_ring_.push_back(x);
//        y_ring_.push_back(0);
//        return 0;
//    }
//}
//
//double fir_filter::filter_opt(double x)
//{
//    x_ring_.push_back(x);
//    double y;
//    switch (nr_coefficients_)
//    {
//    case 5:
//        y = rolling_fir_filter_5_(d, n, x_ring_.front_elements(), y_ring_.front_elements());
//        break;
//    case 4:
//        y = rolling_fir_filter_4_(d, n, x_ring_.front_elements(), y_ring_.front_elements());
//        break;
//    case 3:
//        y = rolling_fir_filter_3_(d, n, x_ring_.front_elements(), y_ring_.front_elements());
//        break;
//    case 2:
//        y = rolling_fir_filter_2_(d, n, x_ring_.front_elements(), y_ring_.front_elements());
//        break;
//    }
//    y_ring_.push_back(y);
//    y_ring_.pop_front();
//    x_ring_.pop_front();
//    return y;
//}
//
//void fir_filter::init_history_values(double x, int nr_samples)
//{
//    for (int i = 0; i < 4 * nr_samples; ++i)
//        filter(x);
//}
//
//i_fir_filter* i_fir_filter::new_instance(const double *n, const double *d, size_t nr_coefficients)
//{
//    return new fir_filter(n, d, nr_coefficients);
//}
//
//void i_fir_filter::delete_instance(i_fir_filter* instance)
//{
//    delete((fir_filter*)instance);
//}

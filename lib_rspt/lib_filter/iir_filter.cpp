/***************************************************************************
* Copyright 2024 Tamas Levente Kis - tamkis@gmail.com                      *
*                                                                          *
* Licensed under the Apache License, Version 2.0 (the "License");          *
* you may not use this file except in compliance with the License.         *
* You may obtain a copy of the License at                                  *
*                                                                          *
*     http://www.apache.org/licenses/LICENSE-2.0                           *
*                                                                          *
* Unless required by applicable law or agreed to in writing, software      *
* distributed under the License is distributed on an "AS IS" BASIS,        *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
* See the License for the specific language governing permissions and      *
* limitations under the License.                                           *
***************************************************************************/
#include <inttypes.h>
#include <string.h>
#include <vector>

#include "../../FifoBuffer.h"
#include "../iir_filter.h"

static inline double rolling_iir_filter_5_(const double* n, const double* d, const double* xz, const double* yz)
{
    return n[0] * xz[4] + n[1] * xz[4 - 1] + n[2] * xz[4 - 2] + n[3] * xz[4 - 3] + n[4] * xz[4 - 4] - d[1] * yz[4 - 1] - d[2] * yz[4 - 2] - d[3] * yz[4 - 3] - d[4] * yz[4 - 4];
}

static inline double rolling_iir_filter_4_(const double* n, const double* d, const double* xz, const double* yz)
{
    return n[0] * xz[3] + n[1] * xz[3 - 1] + n[2] * xz[3 - 2] + n[3] * xz[3 - 3] - d[1] * yz[3 - 1] - d[2] * yz[3 - 2] - d[3] * yz[3 - 3];
}

static inline double rolling_iir_filter_3_(const double* n, const double* d, const double* xz, const double* yz)
{
    return n[0] * xz[2] + n[1] * xz[2 - 1] + n[2] * xz[2 - 2] - d[1] * yz[2 - 1] - d[2] * yz[2 - 2];
}

static inline double rolling_iir_filter_2_(const double* n, const double* d, const double* xz, const double* yz)
{
    return n[0] * xz[1] + n[1] * xz[1 - 1] - d[1] * yz[1 - 1];
}

class iir_filter: public i_iir_filter
{
    fifo_buffer<double> x_ring_;
    fifo_buffer<double> y_ring_;
    double n[5];
    double d[5];
    size_t nr_coefficients_;
public:

    iir_filter(const double *an, const double *ad, size_t nr_coefficients)
        : x_ring_(0),
          y_ring_(0),
          nr_coefficients_(nr_coefficients)
    {
        memcpy(n, an, nr_coefficients * sizeof(double));
        memcpy(d, ad, nr_coefficients * sizeof(double));
    }

    double filter(double x)
    {
        if (x_ring_.size() >= (nr_coefficients_ - 1))
            return filter_opt(x);
        else
        {
            x_ring_.push_back(x);
            y_ring_.push_back(0);
            return 0;
        }
    }

    double filter_opt(double x)
    {
        x_ring_.push_back(x);
        double y;
        switch (nr_coefficients_)
        {
        case 5:
            y = rolling_iir_filter_5_(d, n, x_ring_.front_elements(), y_ring_.front_elements());
            break;
        case 4:
            y = rolling_iir_filter_4_(d, n, x_ring_.front_elements(), y_ring_.front_elements());
            break;
        case 3:
            y = rolling_iir_filter_3_(d, n, x_ring_.front_elements(), y_ring_.front_elements());
            break;
        case 2:
            y = rolling_iir_filter_2_(d, n, x_ring_.front_elements(), y_ring_.front_elements());
            break;
        }
        y_ring_.push_back(y);
        y_ring_.pop_front();
        x_ring_.pop_front();
        return y;
    }

    void init_history_values(double x, int nr_samples)
    {
        for (int i = 0; i < 4 * nr_samples; ++i)
            filter(x);
    }

    virtual ~iir_filter() = default;
};

i_iir_filter* i_iir_filter::new_instance(const double *n, const double *d, size_t nr_coefficients)
{
    return new iir_filter(n, d, nr_coefficients);
}

void i_iir_filter::delete_instance(i_iir_filter* instance)
{
    delete((iir_filter*)instance);
}

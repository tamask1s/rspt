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
#include <iostream>
#include <cmath>

using namespace std;

#include "../filter.h"

static inline void rolling_iir_filter_5_(const double* d, const double* n, const double* xz, double* yz)
{
    yz[0] = d[0] * xz[0] + d[1] * xz[1] + d[2] * xz[2] + d[3] * xz[3] + d[4] * xz[4] - n[1] * yz[1] - n[2] * yz[2] - n[3] * yz[3] - n[4] * yz[4];
}

static inline void rolling_iir_filter_4_(const double* d, const double* n, const double* xz, double* yz)
{
    yz[0] = d[0] * xz[0] + d[1] * xz[1] + d[2] * xz[2] + d[3] * xz[3] - n[1] * yz[1] - n[2] * yz[2] - n[3] * yz[3];
}

static inline void rolling_iir_filter_3_(const double* d, const double* n, const double* xz, double* yz)
{
    yz[0] = d[0] * xz[0] + d[1] * xz[1] + d[2] * xz[2] - n[1] * yz[1] - n[2] * yz[2];
}

static inline void rolling_iir_filter_2_(const double* d, const double* n, const double* xz, double* yz)
{
    yz[0] = d[0] * xz[0] + d[1] * xz[1] - n[1] * yz[1];
}

class iir_filter: public i_filter
{
    std::vector<double> x_ring_;
    std::vector<double> y_ring_;
    double n[5];
    double d[5];
    size_t nr_coefficients_;
public:

    iir_filter(const double *an, const double *ad, size_t nr_coefficients)
        :x_ring_(nr_coefficients, 0.0),
         y_ring_(nr_coefficients, 0.0),
         nr_coefficients_(nr_coefficients)
    {
        memcpy(n, an, nr_coefficients * sizeof(double));
        memcpy(d, ad, nr_coefficients * sizeof(double));
    }

    double filter(double x)
    {
        for (int i = nr_coefficients_ - 1; i > 0; --i)
        {
            x_ring_[i] = x_ring_[i - 1];
            y_ring_[i] = y_ring_[i - 1];
        }
        x_ring_[0] = x;
        y_ring_[0] = d[0] * x_ring_[0];
        for (unsigned int i = 1; i < nr_coefficients_; ++i)
        {
            y_ring_[0] += d[i] * x_ring_[i];
            y_ring_[0] -= n[i] * y_ring_[i];
        }
        return y_ring_[0];
    }

    double filter_opt(double x)
    {
        for (int i = nr_coefficients_ - 1; i > 0; --i)
        {
            x_ring_[i] = x_ring_[i - 1];
            y_ring_[i] = y_ring_[i - 1];
        }
        x_ring_[0] = x;
        switch (nr_coefficients_)
        {
        case 5:
            rolling_iir_filter_5_(d, n, x_ring_.data(), y_ring_.data());
            break;
        case 4:
            rolling_iir_filter_4_(d, n, x_ring_.data(), y_ring_.data());
            break;
        case 3:
            rolling_iir_filter_3_(d, n, x_ring_.data(), y_ring_.data());
            break;
        case 2:
            rolling_iir_filter_2_(d, n, x_ring_.data(), y_ring_.data());
            break;
        default:
            break;
        }
        return y_ring_[0];
    }

    void init_history_values(double x, int nr_samples)
    {
        for (int i = 0; i < 4 * nr_samples; ++i)
            filter(x);
    }

    virtual ~iir_filter() = default;
};

i_filter* i_filter::new_iir(const double *n, const double *d, size_t nr_coefficients)
{
    return new iir_filter(n, d, nr_coefficients);
}

void i_filter::delete_iir(i_filter* instance)
{
    delete((iir_filter*)instance);
}

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

#include "../lib_ring_buffer/ring_buffers.h"
#include "../iir_filter.h"

using namespace std;

class fir_filter: public i_filter
{
    fifo_buffer<double> x_ring_;
    vector<double> kernel_;
    size_t kernel_size_;
public:

    fir_filter(const double *kernel, size_t kernel_size)
        : x_ring_(0),
          kernel_size_(kernel_size)
    {
        kernel_.resize(kernel_size);
        memcpy(kernel_.data(), kernel, kernel_size * sizeof(double));
    }

    double filter(double x)
    {
        if (x_ring_.size() == kernel_size_)
            return filter_opt(x);
        else
        {
            x_ring_.push_back(x);
            return 0;
        }
    }

    double filter_opt(double x)
    {
        x_ring_.push_back(x);
        x_ring_.pop_front();
        double y = 0;
        for (unsigned int i = 0; i < kernel_size_; ++i)
            y += x_ring_[i] * kernel_[i];
        return y;
    }

    void init_history_values(double x, int nr_samples)
    {
        for (unsigned int i = 0; i < kernel_size_; ++i)
            filter(x);
    }

    virtual ~fir_filter() = default;
};

i_filter* i_filter::new_fir(const double *kernel, size_t kernel_size)
{
    return new fir_filter(kernel, kernel_size);
}

void i_filter::delete_fir(i_filter* instance)
{
    delete((fir_filter*)instance);
}

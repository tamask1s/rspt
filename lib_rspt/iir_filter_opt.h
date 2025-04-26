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

struct iir_filter_2nd_order
{
    double xz[3];
    double yz[3];
    vector<double> n;
    vector<double> d;
    static const size_t nr_coefficients_ = 3;
    iir_filter_2nd_order()
        : n(3, 0.0), d(3, 0.0)
    {
        memset(xz, 0, nr_coefficients_ * sizeof(double));
        memset(yz, 0, nr_coefficients_ * sizeof(double));
    }

    inline double filter(double x)
    {
        xz[2] = xz[1];
        yz[2] = yz[1];
        xz[1] = xz[0];
        yz[1] = yz[0];
        xz[0] = x;
        yz[0] = d[0] * xz[0] + d[1] * xz[1] + d[2] * xz[2] - n[1] * yz[1] - n[2] * yz[2];
        return yz[0];
    }

    void init_history_values(double x, int nr_samples)
    {
        for (int i = 0; i < 4 * nr_samples; ++i)
            filter(x);
    }
};

struct iir_filter_1st_order
{
    double xz[2];
    double yz[2];
    vector<double> n;
    vector<double> d;
    static const size_t nr_coefficients_ = 2;
    iir_filter_1st_order()
        : n(2, 0.0), d(2, 0.0)
    {
        memset(xz, 0, nr_coefficients_ * sizeof(double));
        memset(yz, 0, nr_coefficients_ * sizeof(double));
    }

    inline double filter(double x)
    {
        xz[1] = xz[0];
        yz[1] = yz[0];
        xz[0] = x;
        yz[0] = d[0] * xz[0] + d[1] * xz[1] - n[1] * yz[1];
        return yz[0];
    }

    void init_history_values(double x, int nr_samples)
    {
        for (int i = 0; i < 4 * nr_samples; ++i)
            filter(x);
    }
};

struct iir_filter_4th_order
{
    double xz[5];
    double yz[5];
    vector<double> n;
    vector<double> d;
    size_t nr_coefficients_ = 5;

    iir_filter_4th_order()
        : n(5, 0.0), d(5, 0.0)
    {
        memset(xz, 0, nr_coefficients_ * sizeof(double));
        memset(yz, 0, nr_coefficients_ * sizeof(double));
    }

    inline double filter(double x)
    {
        for (int i = nr_coefficients_ - 1; i > 0; --i)
        {
            xz[i] = xz[i - 1];
            yz[i] = yz[i - 1];
        }
        xz[0] = x;
        yz[0] = d[0] * xz[0] + d[1] * xz[1] + d[2] * xz[2] + d[3] * xz[3] + d[4] * xz[4] - n[1] * yz[1] - n[2] * yz[2] - n[3] * yz[3] - n[4] * yz[4];
        return yz[0];
    }

    void init_history_values(double x, int nr_samples)
    {
        for (int i = 0; i < 4 * nr_samples; ++i)
            filter(x);
    }
};

class delay
{
    vector<double> history;
    size_t nr_samples_;
public:
    delay(size_t nr_samples)
        :history(nr_samples, 0.0),
         nr_samples_(nr_samples)
    {}
    inline double get_delayed(const double new_sample)
    {
        double res = history[nr_samples_ - 1];
        for (int i = nr_samples_ - 1; i > 0; --i)
            history[i] = history[i - 1];
        history[0] = new_sample;
        return res;
    }
};

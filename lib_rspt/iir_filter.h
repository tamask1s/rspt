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

/** example coefficients:

double n[] = {1.00000000000, -3.14332095199, 3.70064088865, -1.97083923944, 0.41351972908}; /// BP 0.4-200
double d[] = {0.06722876941, 0.00000000000, -0.13445753881, 0.00000000000, 0.06722876941};

double n[] = {1.00000000000, -1.56101807580, 0.64135153806}; /// LP 100Hz
double d[] = {0.02008336556, 0.04016673113, 0.02008336556};

double n[] = {1.00000000000, -1.99955571171, 0.99955581039}; /// HP 0.1Hz
double d[] = {0.99977788053, -1.99955576105, 0.99977788053};

double n[] = {1.00000000000, -1.99822284729, 0.99822442503}; /// HP 0.4Hz
double d[] = {0.99911181808, -1.99822363616, 0.99911181808};

double n[] = {1.00000000000, -1.99555712435, 0.99556697207}; /// HP 1Hz
double d[] = {0.99778102410, -1.99556204821, 0.99778102410};

double n[] = {1.00000000000, -1.56101807580, 0.64135153806}; /// HP 100Hz
double d[] = {0.80059240346, -1.60118480693, 0.80059240346};

double fir_kernel[] = {0.111, 0.111 * 2, 0.111 * 3, 0.111 * 2, 0.111};  /// some LP
double fir_kernel[] = {-0.2, -0.3, -0.5, 0, 0.5, 0.3, 0.2};  /// some HP

*/

class i_filter
{
public:
    virtual double filter(double x) = 0;
    virtual double filter_opt(double x) = 0;
    virtual void init_history_values(double x, int nr_samples = 2000) = 0;

    static i_filter* new_iir(const double *n, const double *d, size_t nr_coefficients);
    static void delete_iir(i_filter* instance);

    static i_filter* new_fir(const double *kernel, size_t kernel_size);
    static void delete_fir(i_filter* instance);
};

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

class i_signal_packer
{
public:
    virtual void compress(const unsigned char* src, unsigned char* dst, size_t dst_max_len, size_t& dst_len) = 0;
    virtual void decompress(const unsigned char* src, size_t& src_len, unsigned char* dst) = 0;

    static i_signal_packer* new_xdelta_hzr(size_t bytes_per_channel, size_t nr_of_channels, size_t nr_of_samples_in_each_channel);
    static void delete_xdelta_hzr(i_signal_packer* instance);

    static i_signal_packer* new_hzr(size_t bytes_per_channel, size_t nr_of_channels, size_t nr_of_samples_in_each_channel);
    static void delete_hzr(i_signal_packer* instance);

    static i_signal_packer* new_dct(size_t bytes_per_channel, size_t nr_of_channels, size_t nr_of_samples_in_each_channel);
    static void delete_dct(i_signal_packer* instance);

    static i_signal_packer* new_hadamard(size_t bytes_per_channel, size_t nr_of_channels, size_t nr_of_samples_in_each_channel);
    static void delete_hadamard(i_signal_packer* instance);
};

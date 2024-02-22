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

/**
 * @brief Abstract class for signal compression and decompression.
 *
 * This class defines the interface for signal compression and decompression algorithms.
 * The class's compression algorithms operate on fixed-size input data.
 * The efficiency of the packers requires knowledge of the internal structures of the
 * data to be compressed, therefore during object initialization we provide not only
 * the (fixed) datasize, but both BYTESPERSAMPLE, nr_channels and nr_samples needs to
 * be  specified. The total input datasize will be the product of these in bytes,
 * BYTESPERSAMPLE being the number of bytes one sample of one single channel accupies,
 * while nr_samples is the number of samples for each channel.
 */
class i_signal_packer
{
public:
    /**
     * @brief Compresses the input signal.
     *
     * Compresses the input signal stored in 'src' and writes the compressed data to 'dst'.
     * The maximum allowed size for the compressed data is specified by 'dst_max_len'.
     * The actual size of the compressed data is written to 'dst_len'.
     *
     * @param src The input signal data to be compressed.
     * @param dst The buffer where the compressed data will be written.
     * @param dst_max_len The maximum allowed size for the compressed data.
     * @param dst_len The actual size of the compressed data.
     */
    virtual void compress(const unsigned char* src, unsigned char* dst, size_t dst_max_len, size_t& dst_len) = 0;

    /**
     * @brief Decompresses the input signal.
     *
     * Decompresses the input signal stored in 'src' and writes the decompressed data to 'dst'.
     * The size of the compressed data is written in 'src_len' by the called function, and it doesn't
     * needs to be provided by the caller.
     *
     * @param src The compressed signal data to be decompressed.
     * @param src_len The size of the compressed data, returned to the caller.
     * @param dst The buffer where the decompressed data will be written.
     */
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

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

class signal_packer_int32_base: public i_signal_packer
{
protected:
    tensor_i32 enc_;
    tensor_ui8 serialized_;
    /// nr_bytes_to_compress In most of the cases the encoded data should fit in 24 bits. Though, theroetical prove or testing with high-amplitudes should be done to be sure.
    virtual void compress_i32(const unsigned char* src, unsigned char* dst, size_t dst_max_len, size_t& dst_len, uint8_t compression_method, unsigned int nr_bytes_to_compress, const unsigned char* header = 0, size_t header_size = 0);
    virtual void decompress_i32(const unsigned char* src, size_t& src_len, unsigned char* dst, uint8_t& compression_method, unsigned int nr_bytes_to_compress, unsigned char* header = 0, size_t header_size = 0);
};

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

#ifndef RING_BUFFERS
#define RING_BUFFERS

template<typename TType>
class continuous_ring
{
    TType* mData;
    unsigned int mRealSize;
    unsigned int mShift;
protected:
    unsigned int mSize;
public:
    TType* mShiftedData;
    continuous_ring(size_t a_size)
    {
        mRealSize = a_size * 2 + 1;
        mData = new TType[mRealSize];
        mShiftedData = mData;
        mSize = 0;
        mShift = 0;
    }
    ~continuous_ring()
    {
        delete[] mData;
    }
    bool empty() const
    {
        return mSize ? 0 : 1;
    }
    size_t size() const
    {
        return mSize;
    }
    inline void push_back(const TType& aElement)
    {
        push_elements_back(&aElement, 1);
    }
    inline void push_elements_back(const TType* aElements, size_t aNrElements)
    {
        if (mShift + mSize + aNrElements > mRealSize)
        {
            if ((mSize + aNrElements <= mRealSize) && (mShift > (mRealSize >> 1)) && ((mRealSize >> 1) >= aNrElements))
                memmove((void*)mData, mShiftedData, mSize * sizeof(TType));
            else
            {
                mRealSize <<= 1;
                if (mRealSize < aNrElements + mSize)
                    mRealSize += aNrElements;
                TType* lNewData = new TType[mRealSize];
                memcpy((void*)lNewData, mShiftedData, mSize * sizeof(TType));
                delete[] mData;
                mData = lNewData;
            }
            mShift = 0;
            mShiftedData = mData;
        }
        memcpy((void*)(mShiftedData + mSize), aElements, aNrElements * sizeof(TType));
        mSize += aNrElements;
    }
    TType* enlarge_back(size_t aNrElements)
    {
        if (mShift + mSize + aNrElements > mRealSize)
        {
            if ((mSize + aNrElements <= mRealSize) && (mShift > (mRealSize >> 1)) && ((mRealSize >> 1) >= aNrElements))
                memmove(mData, mShiftedData, mSize * sizeof(TType));
            else
            {
                mRealSize <<= 1;
                if (mRealSize < aNrElements + mSize)
                    mRealSize += aNrElements;
                TType* lNewData = new TType[mRealSize];
                memcpy(lNewData, mShiftedData, mSize * sizeof(TType));
                delete[] mData;
                mData = lNewData;
            }
            mShift = 0;
            mShiftedData = mData;
        }
        mSize += aNrElements;
        return mShiftedData + mSize - aNrElements;
    }
    void clear()
    {
        mShiftedData = mData;
        mSize = 0;
        mShift = 0;
    }
    TType& operator[] (size_t aIndex) const
    {
        return mShiftedData[aIndex];
    }
    TType& front() const
    {
        return mShiftedData[0];
    }
    TType* front_elements() const
    {
        return mShiftedData;
    }
    TType& back() const
    {
        return mShiftedData[mSize - 1];
    }
    void pop_front()
    {
        if (mSize)
        {
            ++mShiftedData;
            ++mShift;
            --mSize;
        }
    }
    void pop_back()
    {
        if (mSize)
            --mSize;
    }
    void pop_elements_front(size_t aNrElements)
    {
        if (mSize >= aNrElements)
        {
            mShiftedData += aNrElements;
            mShift += aNrElements;
            mSize -= aNrElements;
        }
    }
    void pop_elements_back(size_t aNrElements)
    {
        if (mSize >= aNrElements)
            mSize -= aNrElements;
    }
};

template<size_t nr_max_packets = 100>
class io_buffer
{
    int packet_bytes_;
    int nr_max_packets_;
    std::vector<uint8_t> buffer_;
    int it_read = 0;
    int it_write = 0;
    int it_write_last = 0;
    volatile uint8_t rw_states_[nr_max_packets];

public:
    io_buffer(int packet_size)
        : packet_bytes_(packet_size),
          nr_max_packets_(nr_max_packets)
    {
        buffer_.resize(packet_bytes_ * nr_max_packets_);
        for (int i = 0; i < nr_max_packets_; ++i)
            rw_states_[i] = 0;
    }

    uint8_t* get_next_filled_address()
    {
        uint8_t* res = 0;
        if (rw_states_[it_read] == 2)
        {
            res = &buffer_[it_read * packet_bytes_];
            rw_states_[it_read] = 3;
            ++it_read;
            if (it_read == nr_max_packets_)
                it_read = 0;
        }
        return res;
    }

    uint8_t* get_next_address_to_fill()
    {
        uint8_t* res = 0;
        if (rw_states_[it_write] == 0 || rw_states_[it_write] == 3)
        {
            if (rw_states_[it_write_last] == 1)
                rw_states_[it_write_last] = 2;
            res = &buffer_[it_write * packet_bytes_];
            rw_states_[it_write] = 1;
            it_write_last = it_write;
            ++it_write;
            if (it_write == nr_max_packets_)
                it_write = 0;
        }
        return res;
    }
};

#endif /// RING_BUFFERS

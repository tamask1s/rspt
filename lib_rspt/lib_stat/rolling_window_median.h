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
 * Rolling median calculator using std::multiset.
 * Maintains a fixed-size sliding window and calculates the median
 * efficiently by moving an iterator rather than recalculating every time.
 * The RollingWindowMedian class computes the median of a rolling window
 * of fixed size over a stream of data. A "rolling window" maintains
 * the last 'size_' elements and discards the oldest one as new elements are added.
 * The median shall stay correct even for duplicate values.
 * Insertion uses only operations with complexities independent of size_,
 * resulting in an overall complexity of O(1). (However, for larger sizes, the memory
 * operations are slightly more costly, resulting in a slight dependency on size_.)
 */
#define ODD (heap_.size() % 2)
template<typename T>
class rolling_window_median
{
private:
    size_t size_;
    std::multiset<T> heap_;
    continuous_ring<typename std::multiset<T>::iterator> iterator_ring_;
    typename std::multiset<T>::iterator median_it_;

public:
    rolling_window_median(size_t size)
        :size_(size),
         iterator_ring_(size)
    {}

    /**
     * Inserts a new value into the sliding window,
     * and removes the oldest value (if necessary).
     * Updates the median iterator.
     * @param value: The new value to be added to the sliding window.
     * @return: The current median after inserting the value.
     */
    T insert(T value)
    {
        iterator_ring_.push_back(heap_.insert(value));

        if (iterator_ring_.size() > size_)
        {
            const auto& oldest_it_ = iterator_ring_.front();
            const T median_value_ = *median_it_; // copy is more optimal for bigger window sizes
            const T oldest_value_ = *oldest_it_;

            if (value < median_value_)
            {
                if (oldest_value_ < median_value_)
                {
                    heap_.erase(oldest_it_);
                }
                else if (ODD)
                {
                    if (oldest_value_ == median_value_)
                        ++median_it_;
                    heap_.erase(oldest_it_);
                    if (!ODD)
                        --median_it_;
                }
                else
                {
                    --median_it_;
                    if (oldest_value_ == *median_it_)
                        ++median_it_;
                    heap_.erase(oldest_it_);
                    if (!ODD)
                        --median_it_;
                }
            }
            else
            {
                if ((value > median_value_) && (oldest_value_ > median_value_))
                {
                    heap_.erase(oldest_it_);
                }
                else if (ODD)
                {
                    ++median_it_;
                    if (oldest_value_ <= *median_it_)
                        ++median_it_;
                    heap_.erase(oldest_it_);
                    if (!ODD)
                        --median_it_;
                }
                else
                {
                    if (oldest_value_ <= median_value_)
                        ++median_it_;
                    heap_.erase(oldest_it_);
                    if (!ODD)
                        --median_it_;
                }
            }
            iterator_ring_.pop_front();
        }
        else if (heap_.size() > 1)
        {
            if (value < *median_it_)
            {
                if (!ODD)
                    --median_it_;
            }
            else
            {
                if (ODD)
                    ++median_it_;
            }
        }
        else
            median_it_ = heap_.begin();

        if (!ODD)
            return (*median_it_ + *next(median_it_)) / 2.0;
        else
            return *median_it_;
    }
};
